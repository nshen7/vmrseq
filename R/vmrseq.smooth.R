#' Title
#'
#' @param SE
#' @param residSmooth
#' @param bpWindow
#' @param meanSmooth
#' @param bpSpan
#' @param minInSpan
#'
#' @return
#' @export
#'
#' @examples

vmrseq.smooth <- function(
    SE,
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

  gr <- granges(SE)
  M <- assays(SE)[[1]]

  values(gr)$meth <- rowSums(M, na.rm = T)
  values(gr)$total <- rowSums(M>=0, na.rm = T)

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
  mean_meth <- NULL; var <- NULL
  chrs <- as.character(unique(seqnames(gr)))
  for (chromosome in chrs) {
    if (verbose) message("...Chromosome ",
                         paste(chromosome, collapse = ", "), ": ",
                         appendLF = FALSE)

    t1 <- proc.time() # time point 1
    gr_chr <- subset(gr, seqnames(gr) == chromosome)
    M_chr <- M[seqnames(gr) == chromosome, ]
    if (length(gr_chr) < minNumCR) {message("No candidates found."); next}

    # Locfit smooth on fractional methylation if meanSmooth==TRUE
    origin_mean_chr <- gr_chr$meth / gr_chr$total
    if (meanSmooth) {
      fit_chr <- smoothMF(x = start(gr_chr), y = origin_mean_chr,
                          chr = chromosome,
                          weights = gr_chr$total, # across-cell coverage as weights
                          minInSpan = minInSpan, bpSpan = bpSpan,
                          verbose = verbose,
                          parallel = parallel)
      mean_meth_chr <- fit_chr$fitted %>% pmax(0) %>% pmin(1)
      # Keep the raw mean_meth if not smoothed
      ind <- which(!fit_chr$is_smooth)
      if (length(ind) > 0) mean_meth_chr[ind] <- origin_mean_chr[ind]
      if (verbose) message("MF smoothed. ", appendLF = FALSE)
    } else {
      mean_meth_chr <- origin_mean_chr
    }
    mean_meth <- c(mean_meth, mean_meth_chr)

    # Compute variance relative to mean_meth
    var_chr <- computeVar(gr = gr_chr,
                          M = M_chr,
                          mean_meth = mean_meth_chr,
                          bpWindow = bpWindow,
                          parallel = parallel)
    var <- c(var, var_chr)

    t2 <- proc.time()
    if (verbose) message("Variance computed (", round((t2 - t1)[3]/60, 2), " min). ")
  } # end looping chromosome

  if (meanSmooth) values(gr)$smoothed_mean <- mean_meth
  values(gr)$var <- var

  return(gr)
}

