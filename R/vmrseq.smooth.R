#' Title
#'
#' @param SE
#' @param sparseNAdrop
#' @param residSmooth
#' @param bpWindow
#' @param meanSmooth
#' @param bpSpan
#' @param minInSpan
#' @param verbose
#' @param BPPARAM
#'
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom stats fitted median
#' @importFrom gamlss.dist dZIBB dBB
#' @importFrom locfit locfit lp
#' @importFrom DelayedArray rowSums colSums rowMeans colMeans is_sparse
#' @importFrom recommenderlab dropNA2matrix
#' @import dplyr
#' @import GenomicRanges
#'
#' @return
#' @export
#'
#' @examples

vmrseq.smooth <- function(
    SE, sparseNAdrop = is_sparse(assays(SE)[[1]]),
    residSmooth = TRUE, bpWindow = 2000, # param for individual-cell methylation residual smoother
    meanSmooth = FALSE, bpSpan = 0, minInSpan = 0, # params for across-cell mean methylation smoother
    verbose = TRUE, BPPARAM = bpparam()
) {

  # TODO: all the other sanity checks
  if (meanSmooth & bpSpan<=0 & minInSpan<=0)
    stop("If mean methylation need to be smoothed, at least one of 'bpSpan' and 'minInSpan' should be positive (integer) number.")

  for (chromosome in unique(seqnames(SE))) {
    SE_chr <- subset(SE, seqnames(SE) == chromosome)
    if (min(diff(start(SE_chr))) < 2)
      stop("There exists at least 2 rows with position difference less than 2 bp.")
  }

  # TODO: report data dimensions

  gr <- granges(SE)
  M <- assays(SE)[[1]]

  if (!sparseNAdrop) {
    values(gr)$meth <- rowSums(M, na.rm = T)
    values(gr)$total <- rowSums(M >= 0, na.rm = T)
  } else {
    values(gr)$meth <- as.integer(round(rowSums(M)))
    values(gr)$total <- rowSums(M > 0)
  }


  if (min(values(gr)$total) < 3)
    warning("We suggest removing CpG sites with across-cell coverage lower than 3 before running vmrseq.")

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

    # Locfit smooth on fractional methylation if meanSmooth==TRUE
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

