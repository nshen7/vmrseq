bumphunt <- function(gr,
                     cutoff = 0.1,
                     maxGap = 1000, minNumRegion = 5,
                     smooth = T,
                     maxGapSmooth = 2500,
                     minInSpan = 10, bpSpan = 1000,
                     verbose = TRUE,
                     parallel = FALSE) {

  # Compute variance for individual sites
  gr$MF <- with(values(gr), meth / total)
  gr$var <- with(values(gr), MF * (1-MF))

  # Apply smoother on each chromosome serially
  chrs <- as.character(unique(seqnames(gr)))
  for (chromosome in chrs) {
    # Subset one chromosome from gr
    gr_chr <- subset(gr, seqnames(gr) == chromosome)
    if (smooth) {

      weights <- gr_chr$total
      # smooth on one chromosome
      smooth_fit <- smoother(x = start(gr_chr), y = gr_chr$var,
                             chr = chromosome,
                             weights = weights,
                             maxGap = maxGap, minNumRegion = minNumRegion,
                             maxGapSmooth = maxGapSmooth,
                             minInSpan = minInSpan, bpSpan = bpSpan,
                             verbose = verbose,
                             parallel = parallel)

    }

  }

}



smoother <- function(x, y, weights, chr,
                     maxGap = 1000, minNumRegion = 5,
                     maxGapSmooth = 2500,
                     minInSpan = 10, bpSpan = 1000,
                     verbose = TRUE,
                     parallel = FALSE) {


  locfitByCluster <- function(ix) {

    yi <- y[ix]
    xi <- x[ix]
    if (!is.null(weights)) {
      weightsi <- matrix(weights[ix], ncol = 1)
    } else {
      weightsi <- NULL
    }
    clusteri <- clusterC[ix]

    if (is.null((yi)))
      stop("y (var) is missing")
    if (is.null(xi))
      stop("x (pos) is missing")
    if (is.null(clusteri))
      stop("cluster is missing")
    if (is.null(weightsi))
      weightsi <- matrix(1, nrow = nrow(yi), ncol = ncol(yi))

    if (length(ix) >= minNumRegion) {
      sdata <- data.frame(posi = xi, yi = yi,
                          weightsi = weightsi)

      # balance minInSpan and bpSpan
      nn <- minInSpan / length(ix)
      fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan),
                    data = sdata, weights = weightsi, family = "gaussian",
                    maxk = 10000)
      yi <- fitted(fit)
      smoothed <- TRUE
    } else {
      yi <- rep(NA, length(yi))
      smoothed <- FALSE
    }
    return(data.frame(fitted = as.vector(yi), smoothed = smoothed))
  }


  if (!is.null(weights) && is.null(dim(weights)))
    weights <- matrix(weights, ncol = 1)

  t1 <- proc.time()
  chr_len <- rep(chr, each = length(x))
  clusterC <- bumphunter::clusterMaker(factor(chr_len, levels=unique(chr_len)),
                                       x,
                                       assumeSorted = T,
                                       maxGap = maxGapSmooth)

  Indexes <- split(seq(along = clusterC), clusterC)

  if (parallel) {
    ret <- do.call(rbind, bplapply(Indexes, function(idx) locfitByCluster(idx)))
  } else {
    ret <- do.call(rbind, lapply(Indexes, function(idx) locfitByCluster(idx)))
  }

  if (verbose) {
    t2 <- proc.time()
    message("Smoothed (",
            round((t2 - t1)[3]/60, 2), " min). ",
            appendLF = FALSE)
  }

  return(ret) # data.frame with columns: 'fitted' (numeric), 'smoothed' (logical)
}

