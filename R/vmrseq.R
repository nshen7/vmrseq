#' @title Main function for detecting and evaluating significance of VMRs in
#' single-cell bisulfite sequencing data.
#'
#' @description Performs a two-step approach that (1) detects candidate regions,
#' and (2) detect variably methylated region(s) inside each candidate region
#' using a hidden Markov model.
#'
#' @param gr GRanges object containing the chromosome position and methylation
#' values. Should contain two element metadata columns that can be extracted
#' using \code{values(gr)}: `meth` and `total`, indicating number of methylated
#' cells and total number of cells.
#' @param minCov integer scalar value that represents minimum across-cell
#' coverage in QC step. Sites with coverage lower than \code{minCov} are
#' removed. Default value is 5.
#' @param cutoff positive scalar value that represents the cutoff value of
#' variance that is used to discover candidate regions. Default value is 0.10.
#' @param maxGap integer value representing maximum number of basepairs in
#' between neighboring CpGs to be included in the same VMR.
#' @param minNumRegion positive integer that represents the minimum number of
#' CpGs to consider for an VMR as well as a candidate region. Default value is
#' 5. Minimum value is 3.
#' @param smooth logical value that indicates whether or not to smooth the CpG
#' level signal when discovering candidate regions. Defaults to TRUE.
#' @param bpSpan a positive integer that represents the length in base pairs of
#' the smoothing span window if \code{smooth} is TRUE. Default value is 1000.
#' @param minInSpan positive integer that represents the minimum number of CpGs
#' in a smoothing span window if \code{smooth} is TRUE. Default value is 10.
#' @param maxGapSmooth integer value representing maximum number of base pairs
#' in between neighboring CpGs to be included in the same cluster when
#' performing smoothing (should generally be larger than \code{maxGap})
#' @param tp transitProbs object that contains estimated transition
#' probabilities. Can be obtained by the 'estimTransitProbs' function. If
#' \code{tp==NULL}, internal transition probabilities in \code{vmrseq} is used.
#' Default is NULL.
#' @param maxNumMerge positive integer that represents the maximum number of
#' CpGs between two VMRs that can be tolerated when merging VMRs in the same
#' candidate region.Default is 1.
#' @param minNumLong positive integer that represents the minimum number of
#' CpGs to consider for a *long* candidate region. For fast computation use.
#' Default is 20. Long regions will be performed a more thorough search of
#' optimized prevalence value. Minimum value is the \code{minNumRegion}.
#' @param control this sets the control parameters of the outer iterations
#' algorithm. The default setting is the \code{vmrseq.control} function.
#' @param verbose logical value that indicates whether progress messages should
#' be printed to stdout. Defaults value is TRUE.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#'
#' @return a list of two \code{GRanges} object that contains the results of the
#' inference.
#'
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom bumphunter clusterMaker getSegments
#' @importFrom stats fitted median
#' @importFrom gamlss.dist dZIBB dBB
#' @importFrom locfit locfit lp
#' @import dplyr
#' @import GenomicRanges
#'
#'
#' @export
#'
#' @examples
#'
#'
vmrseq <- function(gr,
                   minCov = 5,
                   cutoff = 0.1, # param for CR calling
                   maxGap = 1000, minNumRegion = 5, # params for VMR calling
                   smooth = TRUE, maxGapSmooth = 2500, # params for smoother
                   bpSpan = 1000, minInSpan = 10, # params for smoother
                   tp = NULL,
                   maxNumMerge = 1, minNumLong = 10,
                   control = vmrseq.control(),
                   verbose = TRUE, BPPARAM = bpparam()) {

  if (is.null(cutoff) | length(cutoff) != 1 | cutoff <= 0)
    stop("'cutoff' has to be a postive scalar value.")
  if (minNumRegion < 3)
    stop("'minNumRegion' must be at least 3.")
  if (minNumLong < minNumRegion)
    stop("'minNumLong' must be greater or equal to `minNumRegion`.")
  if (class(gr)[1] != "GRanges")
    stop("'gr' must be a GRanges object.")
  if (is.null(gr$total) | all(gr$total != round(gr$total)))
    stop("'gr' must contain an integer column 'total' in element metadata.")
  if (is.null(gr$meth) | all(gr$meth != round(gr$meth)))
    stop("'gr' must contain an integer column 'meth' in element metadata.")
  if (!is.null(tp) & class(tp) != "transitProbs")
    stop("'tp' must be a transitProbs object. Can be estimated by function `estimTransitProbs`.")

  # Check that GRanges object is sorted
  if (is.unsorted(gr)) {
    message("'gr' is not sorted. Sorting 'gr' now.")
    gr <- sort(gr)
  }

  for (chromosome in unique(seqnames(gr))) {
    gr_chr <- subset(gr, seqnames(gr) == chromosome)
    if (min(diff(start(gr_chr))) < 2)
      stop("There exists at least 2 rows with position difference less than 2 bp.")
  }

  # QC: remove low-coverage sites
  if (minCov > 0 & min(gr$total) < minCov) {
    gr <- subset(gr, gr$total >= minCov)
    message("Removed sites with coverage lower than ", minCov)
  }

  # Register the parallel backend
  BiocParallel::register(BPPARAM)
  backend <- paste0("BiocParallel:", class(bpparam())[1])

  if (bpparam()$workers == 1) {
    if (verbose) {
      mes <- "Using a single core (backend: %s)."
      message(sprintf(mes, backend))
    }
    parallel <- FALSE
  } else {
    if (verbose) {
      mes <- paste0("Parallelizing using %s workers/cores ",
                    "(backend: %s).")
      message(sprintf(mes, bpparam()$workers, backend))
    }
    parallel <- TRUE
  }

  message("Detecting candidate regions with smoothed variance larger than ",
          cutoff, "...")
  # Bumphunt candidate regions. Outputs list of index vectors.
  # Each list element is one CR.
  CRI <- callCandidRegion(
    gr = gr,
    cutoff = cutoff,
    maxGap = maxGap, minNumRegion = minNumRegion,
    smooth = smooth,
    maxGapSmooth = maxGapSmooth,
    minInSpan = minInSpan, bpSpan = bpSpan,
    verbose = verbose,
    parallel = parallel
  )

  if (length(CRI) == 0) {
    message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    return(NULL)
  } else {
    message("Finished calling candidate regions (found ", length(CRI),
            " candidate regions in total).")
  }

  message("Detecting VMRs...")
  t1 <- proc.time()
  # Outputs a GRanges objects with VMR ranges and summary information
  VMRI <- searchVMR(
    gr = gr,
    CRI = CRI,
    maxGap = maxGap, minNumRegion = minNumRegion,
    tp = tp,
    maxNumMerge = maxNumMerge,
    minNumLong = minNumLong,
    control = control,
    verbose = verbose,
    parallel = parallel
  )

  t2 <- proc.time()
  if (length(VMRI) == 0) {
    message("No VMR detected.")
    return(NULL)
  } else {
    message("Finished detecting VMRs (took ",
            round((t2 - t1)[3]/60, 2), " min and ",
            nrow(VMRI), " VMRs found in total).")
  }

  vmr.gr <- indexToGranges(gr = gr, index = VMRI, type = "VMR")
  cr.gr <- indexToGranges(gr = gr, index = CRI, type = "CR")

  return(list(VMRs = vmr.gr, CRs = cr.gr))
}



