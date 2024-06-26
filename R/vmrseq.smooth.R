#' @title Smoothing on single-cell bisulfite sequencing data for the purpose of
#' constructing candidate regions.
#'
#' @description \code{vmrseq.smooth} takes a \code{SummarizedExperiment} object
#'  with information of methylation level of individual cells as input, and
#'  perform a kernel smoother to ‘relative’ methylation levels of individual
#'  cells prior to constructing candidate regions. Purpose of the smoothing is
#'  to adjust for uneven coverage biases and borrow information from nearby sites.
#'  See manuscript for detailed description.
#'
#' @param SE \code{SummarizedExperiment} object with one (and only one) assay that
#'  contains *binary* methylation status of CpG sites in individual cells. We
#'  recommend using output by \code{vmrseq::data.pool} (i.e., an NA-dropped
#'  HDF5-based SummarizedExperiment object) to prevent running out of memory.
#' @param bpWindow positive integer that represents the width (in bp) of
#'  smoothing window. Default value is 2000.
#' @param sparseNAdrop logical value that represents whether the NA values are
#'  droppped in the input \code{SE} object. \code{SE} objects output by
#'  \code{vmrseq::data.pool} are NA dropped. See \code{?vmrseq::data.pool}
#'  for details about NA-dropped representation.
#' @param verbose logical value that indicates whether progress messages
#'  should be printed to stdout. Defaults value is TRUE.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#'  backend. The default option is \code{BiocParallel::bpparam()} which will
#'  automatically creates a cluster appropriate for the operating system.
#'
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom stats fitted median
#' @importFrom gamlss.dist dZIBB dBB
#' @importFrom locfit locfit lp
#' @importFrom DelayedArray rowSums colSums rowMeans colMeans is_sparse
#' @importFrom recommenderlab dropNA2matrix
#' @import GenomicRanges
#' @import SummarizedExperiment
#'
#' @return a \code{GRanges} object that contains the result of smoothing.
#'  The object retains genomic coordinates (chr, start, end) of input CpG
#'  sites, in the same order as in the input \code{SE} object. Three
#'  column are added (on top of original metadata columns for the CpG sites in
#'  \code{SE}, if any):
#'  1. meth: methylated cell count of the CpG
#'  2. total: total (non-missing) cell count of the CpG
#'  3. var: variance computed based on individual-cell smoothed relative methylation levels.
#'
#' @seealso \code{\link{data.pool}}, \code{\link{vmrseq.fit}}
#' @export
#' 
#' @examples
#' # load example data
#' toy.se <- HDF5Array::loadHDF5SummarizedExperiment(system.file("extdata", "toy", package = "vmrseq"))
#' 
#' # preprocessing
#' total <- DelayedArray::rowSums(SummarizedExperiment::assays(toy.se)$M_mat > 0)
#' toy.se <- subset(toy.se, total >= 3)
#' 
#' # run vmrseq.smooth
#' toy.gr <- vmrseq.smooth(toy.se)
#' toy.gr
#'
vmrseq.smooth <- function(
    SE,
    bpWindow = 2000, # param for individual-cell methylation residual smoother
    sparseNAdrop = is_sparse(assays(SE)[[1]]),
    verbose = TRUE, BPPARAM = bpparam()
) {

  # Params for across-cell mean methylation smoother (for experimental purpose)
  meanSmooth <- FALSE # turning off meanSmooth works better in practice
  bpSpan <- 0; minInSpan <- 0

  if (meanSmooth & bpSpan<=0 & minInSpan<=0)
    stop("If mean methylation need to be smoothed, at least one of 'bpSpan' and 'minInSpan' should be positive (integer) number.")

  for (chromosome in unique(seqnames(SE))) {
    SE_chr <- subset(SE, seqnames(SE) == chromosome)
    if (min(diff(start(SE_chr))) < 2)
      stop("There exists at least 2 rows with position difference less than 2 bp.")
  }


  # TODO: remove sites with total_read = 0
  # TODO: report data dimensions

  gr <- granges(SE)
  M <- assays(SE)[[1]]

  if (!sparseNAdrop) {
    values(gr)$meth <- rowSums(M, na.rm = TRUE)
    values(gr)$total <- rowSums(M >= 0, na.rm = TRUE)
  } else {
    values(gr)$meth <- as.integer(round(rowSums(M)))
    values(gr)$total <- rowSums(M > 0)
  }

#
#   if (min(values(gr)$total) < 3)
#     warning("We suggest removing CpG sites with across-cell coverage lower than 3 before running vmrseq.")

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

  # Smoothing starts
  message("Smoothing in progress...")

  # Apply smoother and compute variance on each chromosome serially
  meanMeth <- NULL; var <- NULL
  chrs <- as.character(unique(seqnames(gr)))
  for (chromosome in chrs) {
    if (verbose) message("...Chromosome ",
                         paste(chromosome, collapse = ", "), ": ",
                         appendLF = FALSE)

    t1 <- proc.time() # time point 1
    gr_chr <- subset(gr, seqnames(gr) == chromosome)
    M_chr <- M[seqnames(gr) == chromosome, ]

    # smooth on fractional methylation if meanSmooth==TRUE
    origin_mean_chr <- gr_chr$meth / gr_chr$total
    if (meanSmooth) {
      fit_chr <- smoothMF(x = start(gr_chr), y = origin_mean_chr,
                          chr = chromosome,
                          weights = gr_chr$total, # across-cell coverage as weights
                          minInSpan = minInSpan, bpSpan = bpSpan,
                          verbose = verbose,
                          parallel = parallel)
      meanMeth_chr <- fit_chr %>% pmax(0) %>% pmin(1)
      if (verbose) message("Mean methylation smoothed. ", appendLF = FALSE)
    } else {
      meanMeth_chr <- origin_mean_chr
    }

    # Compute variance relative to meanMeth
    var_chr <- computeVar(gr = gr_chr,
                          M = M_chr,
                          meanMeth = meanMeth_chr,
                          bpWindow = bpWindow,
                          sparseNAdrop = sparseNAdrop,
                          parallel = parallel)
    t2 <- proc.time()
    if (verbose) message("Variance computed (", round((t2 - t1)[3]/60, 2), " min). ")

    meanMeth <- c(meanMeth, meanMeth_chr)
    var <- c(var, var_chr)
  } # end looping chromosome

  if (meanSmooth) values(gr)$smoothed_mean <- meanMeth
  values(gr)$var <- var

  return(gr)
}

