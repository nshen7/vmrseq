#' @title Main function for detecting and evaluating significance of DMRs.
#'
#' @description Performs a two-step approach that (1) detects candidate regions,
#' and (2) detect variably methylated region(s) inside each candidate region
#' using a hidden Markov model.
#'
#' @param gr GRanges object containing the chromosome position and methylation
#' values. Should contain two element metadata columns that can be extracted using
#' `values(gr)`: `meth` and `total`, indicating number of methylated cells and
#' total number of cells.
#' @param cutoff positive scalar value that represents the cutoff value of
#' variance that is used to discover candidate regions. Default value is 0.10.
#' @param minCov integer scalar value that represents minimum across-cell coverage
#' in QC step. Sites with coverage lower than `minCov` are removed. Default value
#' is 5.
#' @param maxGap integer value representing maximum number of basepairs in
#' between neighboring CpGs to be included in the same VMR.
#' @param minNumRegion positive integer that represents the minimum number of
#' CpGs to consider for an VMR as well as a candidate region. Default value is
#' 5. Minimum value is 3.
#' @param minNumLong positive integer that represents the minimum number of
#' CpGs to consider for a *long* candidate region. For fast computation use.
#' Default is 20. Long regions will be performed a more thorough search of
#' optimized prevalence value. Minimum value is the `minNumRegion`.
#' @param transitProbs
#' @param smooth
#' @param bpSpan
#' @param minInSpan
#' @param maxGapSmooth
#' @param verbose
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#'
#' @return a \code{GRanges} object that contains the results of the inference.
#'
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom bumphunter clusterMaker getSegments
#' @import GenomicRanges
#' @import dplyr
#'
#'
#' @export
#'
#' @examples
#'
#'
vmrseq <- function(gr,
                   cutoff = 0.1, # params for CR calling
                   minCov = 5,
                   maxGap = 1000, minNumRegion = 5, # params for VMR calling
                   transitProbs = NULL,
                   minNumLong = 20,
                   smooth = TRUE,
                   maxGapSmooth = 2500,
                   bpSpan = 1000, minInSpan = 10, # params for smoother
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


  # TODO: formality check for 'transitProbs'

  # Check that GRanges object is sorted
  if (is.unsorted(gr)) {
    message("'gr' is not sorted. Sorting 'gr' now.")
    gr <- sort(gr)
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

  # Bump hunting candidate regions. Outputs list of index vectors,
  # each list element is one CR.
  message("Detecting candidate regions with smoothed variance larger than ", cutoff)
  CRI <- callCandiRegions(gr = gr,
                          cutoff = cutoff,
                          maxGap = maxGap, minNumRegion = minNumRegion,
                          smooth = smooth,
                          maxGapSmooth = maxGapSmooth,
                          minInSpan = minInSpan, bpSpan = bpSpan,
                          verbose = verbose,
                          parallel = parallel)

  if (length(CRI) == 0) {
    message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    return(NULL)
  } else {
    message("Finished calling candidate regions. (", length(CRI), " candidate regions found in total)")
  }

  if (parallel) {

  } else {

  }




}
