#' @title Main function for detecting and evaluating significance of DMRs.
#'
#' @description Performs a two-step approach that (1) detects candidate regions,
#' and (2) detect variably methylated region(s) inside each candidate region
#' using a hidden Markov model.
#'
#' @param gr GRanges object containing the chromosome position and methylation
#' values. Should contain two metadata columns: `meth` and `total`, indicating
#' number of methylated cells and total number of cells.
#' @param cutoff positive scalar value that represents the cutoff value of
#' variance that is used to discover candidate regions. Default value is 0.10.
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
#' @param maxGap
#' @param verbose
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#'
#' @return a \code{GRanges} object that contains the results of the inference.
#'
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom bumphunter clusterMaker getSegments
#'
#' @import GenomicRanges
#' @import tidyverse
#'
#' @export
#'
#' @examples
#'
#'
vmrseq <- function(gr, cutoff = 0.1,
                   minNumRegion = 5, minNumLong = 20,
                   transitProbs = NULL,
                   smooth = TRUE, bpSpan = 1000,
                   minInSpan = 10, maxGapSmooth = 2500, maxGap = 1000,
                   verbose = TRUE, BPPARAM = bpparam()) {

  if (is.null(cutoff) | length(cutoff) != 1 | cutoff <= 0)
    stop("'cutoff' has to be a postive scalar value.")
  if (minNumRegion < 3)
    stop("'minNumRegion' must be at least 3.")
  if (minNumLong < minNumRegion)
    stop("'minNumLong' must be greater or equal to `minNumRegion`.")

  # TODO: formality check for 'transitProbs'

  # Check that GRanges object is sorted
  if (is.unsorted(gr)) {
    message("'gr' is not sorted. Sorting 'gr' now.")
    gr <- sort(gr)
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
  # each list element is one CR
  message("Detecting candidate regions with smoothed variance larger than ", cutoff)
  CR_inds <- bumphunt(gr = gr, minInSpan = minInSpan,
                      minNumRegion = minNumRegion, cutoff = cutoff,
                      maxGap = maxGap,
                      smooth = smooth,
                      maxGapSmooth = maxGapSmooth, bpSpan = bpSpan,
                      verbose = verbose,
                      parallel = parallel)

  # check that at least one candidate region was found; if there were none
  # there is no need to go on to VMR detection

  if (length(CRS) == 0) {
    message("No candidate regions pass the cutoff of ", unique(abs(cutoff)))
    return(NULL)
  } else {

  }




}
