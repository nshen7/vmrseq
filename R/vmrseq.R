#' @title Main function for detecting and evaluating significance of DMRs.
#'
#' @description Performs a two-step approach that (1) detects candidate regions,
#' and (2) detect variably methylated region(s) inside each candidate region
#' using a hidden Markov model.
#'
#' @param gr GRanges object containing the chromosome position and methylation
#' values. Should contain two element metadata columns that can be extracted
#' using `values(gr)`: `meth` and `total`, indicating number of methylated cells
#' and total number of cells.
#' @param minCov integer scalar value that represents minimum across-cell
#' coverage in QC step. Sites with coverage lower than `minCov` are removed.
#' Default value is 5.
#' @param cutoff positive scalar value that represents the cutoff value of
#' variance that is used to discover candidate regions. Default value is 0.10.
#' @param maxGap integer value representing maximum number of basepairs in
#' between neighboring CpGs to be included in the same VMR.
#' @param minNumRegion positive integer that represents the minimum number of
#' CpGs to consider for an VMR as well as a candidate region. Default value is
#' 5. Minimum value is 3.
#' @param smooth
#' @param bpSpan
#' @param minInSpan
#' @param maxGapSmooth
#' @param tp transitProbs object that contains estimated transition
#' probabilities. Can be obtained by the 'estimTransitProbs' function. If
#' `tp==NULL`, internal transition probabilities in `vmrseq` is used. Default is
#' NULL.
#' @param maxNumMerge positive integer that represents the maximum number of
#' CpGs between two VMRs that can be tolerated when merging VMRs in the same
#' candidate region.Default is 1.
#' @param minNumLong positive integer that represents the minimum number of
#' CpGs to consider for a *long* candidate region. For fast computation use.
#' Default is 20. Long regions will be performed a more thorough search of
#' optimized prevalence value. Minimum value is the `minNumRegion`.
#' @param control
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
                   minCov = 5,
                   cutoff = 0.1, # param for CR calling
                   maxGap = 1000, minNumRegion = 5, # params for VMR calling
                   smooth = TRUE, maxGapSmooth = 2500, # params for smoother
                   bpSpan = 1000, minInSpan = 10, # params for smoother
                   tp = NULL,
                   maxNumMerge = 1, minNumLong = 20,
                   control = optimize.control(),
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
  CRI <- callCandiRegions(
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
    message("Finished calling candidate regions. (", length(CRI),
            " candidate regions found in total)")
  }

  message("Detecting VMRs...")
  # Outputs a GRanges objects with VMR ranges and summary information
  VMR <- detectVMRs(
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

  if (length(VMR) == 0) {
    message("No VMR detected.")
    return(NULL)
  } else {
    message("Finished calling VMRs. (", length(VMR), " VMRs found in total)")
  }


  return(VMR)
}




#' Auxiliary function as user interface for vmrseq optimization. Typically only
#' used when calling vmrseq function with the option `control`.
#'
#' @param inits vector of numeric values between 0 and 1 representing initial
#' values of \pi_1 shall be taken in optimization algorithm.
#' @param epsilon numeric value representing the convergence upper bound for
#' the algorithm.
#' @param backtrack logical value indicating whether to use backtracking line
#' search to automatically adjust learning rate. Default is TRUE.
#' @param eta a numeric value representing the learning rate in optimization.
#' Default is `ifelse(backtrack, 0.05, 0.005)`.
#' @param maxIter positive integer value representing the maximum number of
#' iterations in optimization algorithm.
#' @return
#' @export
#'
optimize.control <- function(
    inits = c(.2, .5, .8),
    epsilon = 1e-4,
    backtrack = T,
    eta = ifelse(backtrack, 0.05, 0.005),
    maxIter = 100
) {
  return(list(inits = inits, epsilon = epsilon,
              backtrack = backtrack,
              eta = eta, maxIter = maxIter))
}
