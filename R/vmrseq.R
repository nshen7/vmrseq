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
#' @param penalty
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
#' candidate region. Default is 0.
#' @param minNumLong positive integer that represents the minimum number of
#' CpGs to consider for a *long* candidate region. For fast computation use.
#' Default is 20. Long regions will be performed a more thorough search of
#' optimized prevalence value. Minimum value is the \code{minNumRegion}.
#' @param gradient should gradient descent be applied to optimize pi? If
#' set to FALSE, \code{vmrseq.control()$inits} will be used as candidate values
#' for pi, but no gradient descent will be performed.
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
                   minCov = 3,
                   cutoff = 0.1, # param for CR calling
                   penalty = 4, # params for VMR calling
                   maxGap = 1000, minNumCR = 5, minNumVMR = 5, # params for VMR calling
                   smooth = TRUE, maxGapSmooth = 2500, # params for smoother
                   bpSpan = 10*median(diff(start(gr))), minInSpan = 10, # params for smoother
                   tp = NULL,
                   maxNumMerge = 0, minNumLong = 0,
                   gradient = TRUE,
                   control = vmrseq.control(),
                   verbose = TRUE, BPPARAM = bpparam()) {

  if (is.null(cutoff) | length(cutoff) != 1 | cutoff <= 0)
    stop("'cutoff' has to be a postive scalar value.")
  if (length(penalty)!=1 | penalty < 0)
    stop("'penalty' has to be a non-negative scalar value.")
  # if (minNumRegion < 3)
  #   stop("'minNumRegion' must be at least 3.")
  # if (minNumLong < minNumRegion)
  #   stop("'minNumLong' must be greater or equal to `minNumRegion`.")
  if (!is.logical(gradient))
    stop("'gradient' must be a logical value (TRUE or FALSE).")
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
  if (minCov > 0) {
    pct_rm <- round(sum(gr$total<minCov)/length(gr)*100, 1)
    if (pct_rm >= 10 & pct_rm < 30) {
      message("WARNING:
  Consider lowering 'minCov' value since ", pct_rm, "% sites will be removed due to QC.")
      warning("Consider lowering 'minCov' value since ", pct_rm, "% sites are removed due to QC.")
    } else if (pct_rm >= 30) {
      stop("'minCov' value should be lowered since ", pct_rm, "% sites will be removed due to QC.")
    }
    gr <- subset(gr, gr$total >= minCov)
    message("QC: Removed ", pct_rm, "% of sites with coverage lower than ", minCov, ".")
  }

  # Register the parallel backend
  BiocParallel::register(BPPARAM)
  backend <- paste0("BiocParallel:", class(bpparam())[1])

  if (bpparam()$workers == 1) {
    if (verbose) {
      mes <- "Parallel: Using a single core (backend: %s)."
      message(sprintf(mes, backend))
    }
    parallel <- FALSE
  } else {
    if (verbose) {
      mes <- paste0("Parallel: Parallelizing using %s workers/cores ",
                    "(backend: %s).")
      message(sprintf(mes, bpparam()$workers, backend))
    }
    parallel <- TRUE
  }

  message(
    "Step 1: Detecting candidate regions with (smoothed) variance larger than ",
    cutoff, "..."
  )
  # Bumphunt candidate regions. Outputs list of index vectors.
  # Each list element is one CR.
  res_cr <- callCandidRegion(
    gr = gr,
    cutoff = cutoff,
    maxGap = maxGap, minNumCR = minNumCR,
    smooth = smooth,
    maxGapSmooth = maxGapSmooth,
    minInSpan = minInSpan, bpSpan = bpSpan,
    maxNumMerge = maxNumMerge,
    verbose = verbose,
    parallel = parallel
  )

  # Indexes of candidate regions
  CRI <- res_cr$CRI

  cr_index <- rep(NA, length(gr))
  cr_index[unlist(CRI)] <- rep.int(1:length(CRI), lengths(CRI))

  # Add summary stats (smoothed var and CR index) into output
  values(gr) <- cbind(values(gr), res_cr$smooth_fit, cr_index)

  # Percentage of sites in CRs
  pct_incr <- round(sum(lengths(CRI))/length(gr)*100, 2)
  if (is.null(CRI)) {
    message("...No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    return(NULL)
  } else {
    message("...Finished calling candidate regions - found ", length(CRI),
            " candidate regions in total.
  ...", pct_incr,
            "% QC-passed sites are called to be in candidate regions.")
  }

  if (pct_incr <= 10) {
    message("WARNING:
  Consider lowering 'cutoff' since only ", pct_incr,
            "% QC-passed sites are called to be in candidate region.
    ...If not, might induce low statistical power.")
    warning("Consider lowering 'cutoff' since only ", pct_incr,
            "% QC-passed sites are called to be in candidate region.")
  }

  message(
    "Step 2: Detecting VMRs",
    ifelse(penalty > 0, yes = " with", no = " without"),
    " penalty",
    ifelse(penalty > 0, yes = paste0("=", round(penalty, 2)), no = ""),
    "..."
  )
  t1 <- proc.time()

  # Outputs a GRanges objects with VMR ranges and summary information
  VMRI <- searchVMR(
    gr = gr,
    CRI = CRI,
    penalty = penalty,
    maxGap = maxGap, minNumVMR = minNumVMR,
    tp = tp,
    maxNumMerge = maxNumMerge,
    minNumLong = minNumLong,
    gradient = gradient,
    control = control,
    verbose = verbose,
    parallel = parallel
  )

  t2 <- proc.time()

  if (length(VMRI) == 0) {
    message("No VMR detected.")
    return(NULL)
  } else {
    message("...Finished detecting VMRs - took ",
            round((t2 - t1)[3]/60, 2), " min and ",
            nrow(VMRI), " VMRs found in total.
  ...", round(sum(VMRI$end_ind-VMRI$start_ind+1) / length(gr) * 100, 2),
            "% QC-passed sites are called to be in VMRs.")
  }

  vmr.gr <- indexToGranges(gr = gr, index = VMRI, type = "VMR")
  cr.gr <- indexToGranges(gr = gr, index = CRI, type = "CR")

  # # Add summary stats into output
  # hits_vmr <- findOverlaps(gr, vmr.gr) %>% as.data.frame()
  #
  # vmr_index <- rep(NA, length(gr))
  # vmr_index[hits_vmr$queryHits] <- hits_vmr$subjectHits
  # loglik_diff <- rep(NA, length(gr))
  # loglik_diff[hits_vmr$queryHits] <- vmr.gr$loglik_diff[hits_vmr$subjectHits]
  # values(gr) <- cbind(values(gr), vmr_index, loglik_diff)

  VMRI2 <- lapply(1:nrow(VMRI), function(i) VMRI$start_ind[i]:VMRI$end_ind[i]) # list of indices
  vmr_index <- rep(NA, length(gr))
  vmr_index[unlist(VMRI2)] <- rep.int(1:length(VMRI2), lengths(VMRI2))
  values(gr) <- cbind(values(gr), vmr_index, values(vmr.gr)[vmr_index,])

  return(list(gr_qced = gr, VMRs = vmr.gr, CRs = cr.gr))
}



