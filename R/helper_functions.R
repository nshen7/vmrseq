callCandidRegion <- function(gr,
                            cutoff = 0.1,
                            maxGap = 1000, minNumRegion = 5,
                            smooth = T,
                            maxGapSmooth = 2500,
                            minInSpan = 10, bpSpan = 1000,
                            verbose = TRUE,
                            parallel = FALSE) {

  # Compute variance for individual sites
  gr$MF <- gr$meth / gr$total
  gr$var_raw <- gr$MF * (1-gr$MF)

  # Apply smoother on each chromosome serially
  chrs <- as.character(unique(seqnames(gr)))

  fit <- NULL
  for (chromosome in chrs) {

    if (verbose) message("...Chromosome ", paste(chromosome, collapse = ", "), ": ", appendLF = FALSE)

    # Subset one chromosome from gr
    gr_chr <- subset(gr, seqnames(gr) == chromosome)

    # Skip chromosomes that have fewer than minNumRegion loci
    if (length(gr_chr) < minNumRegion){
      message("No candidates found.")
      next
    }

    # Locfit smooth on var
    if (smooth) {
      weights <- gr_chr$total
      # smooth on one chromosome
      fit_chr <- smoother(x = start(gr_chr), y = gr_chr$var_raw,
                          chr = chromosome,
                          weights = weights,
                          maxGap = maxGap, minNumRegion = minNumRegion,
                          maxGapSmooth = maxGapSmooth,
                          minInSpan = minInSpan, bpSpan = bpSpan,
                          verbose = verbose,
                          parallel = parallel)
    } else {
      fit_chr <- data.frame(fitted = rep(NA, length(gr_chr)), smoothed = FALSE)
    }

    # Keep the raw var if not smoothed
    ind <- which(!fit_chr$smoothed)
    if (length(ind) > 0) {
      fit_chr$fitted[ind] <- gr_chr$var_raw[ind]
    } else {
      fit_chr$fitted <- gr_chr$var_raw
    }

    # Concatenate fit_chr from all chromosomes
    fit <- rbind(fit, fit_chr)
  }

  if (nrow(fit) != length(gr)) {
    stop("nrow(fit) does not match with length(gr).")
  } else {
    gr$var <- fit$fitted
  }

  # Call candidate regions
  cluster <- bumphunter::clusterMaker(chr = seqnames(gr),
                                      pos = start(gr),
                                      maxGap = maxGap,
                                      assumeSorted = TRUE)
  Indexes <- bumphunter::getSegments(x = gr$var, f = cluster,
                                     cutoff = cutoff,
                                     assumeSorted = TRUE,
                                     verbose = FALSE)

  upIndex <- Indexes$upIndex
  CRI <- upIndex[lengths(upIndex) >= minNumRegion]
  return(CRI)
} # end of function `callCandidRegion`



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
      stop("y (var_raw) is missing")
    if (is.null(xi))
      stop("x (pos) is missing")
    if (is.null(clusteri))
      stop("cluster is missing")
    if (is.null(weightsi))
      weightsi <- matrix(1, nrow = nrow(yi), ncol = ncol(yi))

    if (length(ix) >= minNumRegion) {
      sdata <- data.frame(posi = xi, yi = yi, weightsi = weightsi)

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
  } # end of function `locfitByCluster`


  if (!is.null(weights) && is.null(dim(weights)))
    weights <- matrix(weights, ncol = 1)

  t1 <- proc.time()
  chr_len <- rep(chr, each = length(x))
  clusterC <- bumphunter::clusterMaker(factor(chr_len, levels=unique(chr_len)),
                                       x,
                                       assumeSorted = TRUE,
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

} # end of function `smoother`



searchVMR <- function(gr,
                      CRI,
                      maxGap = 1000, minNumRegion = 5,
                      tp = NULL,
                      maxNumMerge = 1,
                      minNumLong = 20,
                      control = optimize.control(),
                      verbose = TRUE,
                      parallel = FALSE) {

  # If no `tp` provided, use internal `tp0`
  if (is.null(tp)) tp <- vmrseq:::tp0

  t1 <- proc.time()

  # Pre-computation to ease computational burden
  max_cov <- max(gr$total)
  med_cov <- median(gr$total)
  REFARRAY <- .calRefArray(max_cov = max_cov)
  CHOICEARRAY <- .calChoiceArray(REFARRAY)
  list <- .calMethArray(par_u = .priorParams(med_cov = med_cov, type = "u"),
                        par_m = .priorParams(med_cov = med_cov, type = "m"),
                        max_cov = max_cov)
  METHARRAY <- list$METHARRAY; UNMETHARRAY <- list$UNMETHARRAY

  searchVMRbyRegion <- function(i) {

    ix <- CRI[[i]]
    totals <- gr$total[ix]
    meths <- gr$meth[ix]
    pos <- start(gr[ix])

    res_1g <- .optim1Grp(pos = pos, totals = totals, meths = meths,
                         tp = tp,
                         METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)

    if (length(ix) >= minNumLong) { # long region
      res_2g <- .optim2Grp(pos = pos, totals = totals, meths = meths,
                           tp = tp,
                           inits = control$inits, epsilon = control$epsilon,
                           backtrack = control$backtrack,
                           eta = control$eta, max_iter = control$maxIter,
                           CHOICEARRAY = CHOICEARRAY,
                           METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
    } else { # short region
      res_2g <- .optim2Grp(pos = pos, totals = totals, meths = meths,
                           tp = tp,
                           inits = mean(meths/totals), epsilon = control$epsilon,
                           backtrack = control$backtrack,
                           eta = control$eta, max_iter = control$maxIter,
                           CHOICEARRAY = CHOICEARRAY,
                           METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
    }

    if (res_2g$loglik > res_1g$loglik) {
      vmr_inds <- .callVMR(
        state_seq = res_2g$vit_path[, 1:2],
        min_n = minNumRegion,
        max_n_merge = maxNumMerge
      )
      if (is.null(vmr_inds)) return(NULL)
      else return(data.frame(vmr_inds + ix[1] - 1,
                             cr_name = i,
                             optim_pi = res_2g$optim_pi_1))
    } else {
      return(NULL)
    }
  } # end of function `searchVMRbyRegion`

  if (parallel) {
    VMRI <- do.call(rbind, bplapply(1:length(CRI), function(i) searchVMRbyRegion(i)))
  } else {
    VMRI <- do.call(rbind, lapply(1:length(CRI), function(i) searchVMRbyRegion(i)))
  }

  if (verbose) {
    t2 <- proc.time()
    message("Finished VMR detection (",
            round((t2 - t1)[3]/60, 2), " min). ",
            appendLF = FALSE)
  }

  # return  a data.frame with columns:
  # 'start_ind' (start index of VMR), 'end_ind' (end index of VMR), 'cr_name' (name of CR)
  return(VMRI)

} # end of function `searchVMR`




