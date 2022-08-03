callCandidRegion <- function(gr,
                             cutoff = 0.05,
                             maxGap = 1000, minNumRegion = 5,
                             smooth = T,
                             maxGapSmooth = 2500,
                             minInSpan = 10, bpSpan = 10*median(diff(start(gr))),
                             maxNumMerge = 1,
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
      fit_chr <- data.frame(fitted_var = rep(NA, length(gr_chr)), is_smooth = FALSE)
    }

    # Keep the raw var if not smoothed
    ind <- which(!fit_chr$is_smooth)
    if (length(ind) > 0) fit_chr$fitted_var[ind] <- gr_chr$var_raw[ind]

    # Concatenate fit_chr from all chromosomes
    fit <- rbind(fit, fit_chr)
  }

  # if (nrow(fit) != length(gr)) stop("nrow(fit) does not match with length(gr).")


  # Call candidate regions
  cluster <- bumphunter::clusterMaker(chr = seqnames(gr),
                                      pos = start(gr),
                                      maxGap = maxGap,
                                      assumeSorted = TRUE)
  Indexes <- bumphunter::getSegments(x = fit$fitted_var, f = cluster,
                                     cutoff = cutoff,
                                     assumeSorted = TRUE,
                                     verbose = FALSE)

  upIndex <- Indexes$upIndex
  if (length(upIndex) == 0) upIndex <- NULL

  # Merge CRs if they are closer than `maxNumMerge` bp
  if (!is.null(upIndex)) {
    if (maxNumMerge > 0) {
      for (i in 1:(length(upIndex)-1)) {
        fr <- upIndex[[i]]
        bh <- upIndex[[i+1]]
        if (bh[1] - fr[length(fr)] <= maxNumMerge + 1) {
          combined <- (fr[1]):(bh[length(bh)])
          upIndex[[i]] <- NA
          upIndex[[i+1]] <- combined
        }
      }
      upIndex <- upIndex[!is.na(upIndex)]
    }
    # Only keep candidate regions with more than `minNumRegion` CpGs
    CRI <- upIndex[lengths(upIndex) >= minNumRegion]
    return(list(CRI = CRI, smooth_fit = fit))
  } else {
    return(list(CRI = NULL, smooth_fit = fit))
  }

} # end of function `callCandidRegion`



smoother <- function(x, y, weights, chr,
                     maxGap, minNumRegion,
                     maxGapSmooth,
                     minInSpan, bpSpan,
                     verbose = TRUE,
                     parallel = FALSE) {


  locfitByCluster <- function(idx) {

    yi <- y[idx]
    xi <- x[idx]
    if (!is.null(weights)) {
      weightsi <- weights[idx]
    } else {
      weightsi <- NULL
    }
    clusteri <- clusterC[idx]

    if (is.null((yi)))
      stop("y (var_raw) is missing")
    if (is.null(xi))
      stop("x (pos) is missing")
    if (is.null(clusteri))
      stop("cluster is missing")

    if (length(idx) >= minNumRegion) {
      df <- data.frame(posi = xi, yi = yi, weightsi = weightsi)

      # balance minInSpan and bpSpan
      nn <- minInSpan / length(idx)
      fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan),
                    data = df, weights = weightsi, family = "gaussian",
                    maxk = 10000)
      yi <- fitted(fit)
      is_smooth <- TRUE
    } else {
      yi <- rep(NA, length(yi))
      is_smooth <- FALSE
    }
    return(data.frame(fitted_var = as.vector(yi), is_smooth = is_smooth))
  } # end of function `locfitByCluster`


  # if (!is.null(weights) && is.null(dim(weights)))
  #   weights <- matrix(weights, ncol = 1)

  t1 <- proc.time()
  chr_len <- rep(chr, each = length(x))
  clusterC <- bumphunter::clusterMaker(factor(chr_len, levels=unique(chr_len)),
                                       x,
                                       assumeSorted = TRUE,
                                       maxGap = maxGapSmooth)

  Indexes <- split(seq(along = clusterC), clusterC)

  if (parallel) {
    ret <- do.call(rbind, bplapply(Indexes, function(idx) locfitByCluster(idx)))
    # ret <- do.call(rbind, parallel::mclapply(Indexes, function(idx) locfitByCluster(idx), mc.cores = 6))
  } else {
    ret <- do.call(rbind, lapply(Indexes, function(idx) locfitByCluster(idx)))
  }

  if (verbose) {
    t2 <- proc.time()
    message("Smoothed (", round((t2 - t1)[3]/60, 2), " min). ")
    # message("Smoothed (",
    #         round((t2 - t1)[3]/60, 2), " min). ",
    #         appendLF = FALSE)
  }

  return(ret) # data.frame with columns: 'fitted_var' (numeric), 'is_smooth' (logical)

} # end of function `smoother`



searchVMR <- function(gr,
                      CRI,
                      penalty = penalty,
                      maxGap = 1000, minNumRegion = 5,
                      tp = NULL,
                      maxNumMerge = 1,
                      minNumLong = 10,
                      control = vmrseq.control(),
                      verbose = TRUE,
                      parallel = FALSE) {

  # If no `tp` provided, use internal `tp0`
  if (is.null(tp)) tp <- tp0

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

    if (res_2g$loglik > res_1g$loglik + penalty) {
      vmr_inds <- .callVMR(
        state_seq_2g = res_2g$vit_path[, 1:2],
        min_n = minNumRegion,
        max_n_merge = maxNumMerge
      )
      if (is.null(vmr_inds)) return(NULL)
      else return(data.frame(vmr_inds + ix[1] - 1,
                             cr_index = i,
                             vmr_num_cpg = vmr_inds$end_ind - vmr_inds$start_ind + 1,
                             optim_pi = res_2g$optim_pi_1,
                             n_iter = res_2g$n_iter,
                             loglik_diff = res_2g$loglik - res_1g$loglik - penalty))
    } else {
      return(NULL)
    }
  } # end of function `searchVMRbyRegion`

  if (parallel) {
    VMRI <- do.call(
      rbind,
      bplapply(1:length(CRI), function(i) searchVMRbyRegion(i))
    )
  } else {
    VMRI <- do.call(
      rbind,
      lapply(1:length(CRI), function(i) searchVMRbyRegion(i))
    )
  }

  # return  a data.frame with columns:
  # 'start_ind' (start index of VMR), 'end_ind' (end index of VMR), 'cr_index' (name of CR)
  return(VMRI)

} # end of function `searchVMR`



indexToGranges <- function(gr, index, type) {

  if (!type %in% c("VMR", "CR")) stop("'type' has to be 'VMR' or 'CR'.")

  if (type == "VMR") {
    vmr.gr <- GRanges(
      seqnames = seqnames(gr)[index$start_ind],
      ranges = IRanges(start = start(gr)[index$start_ind],
                       end = end(gr)[index$end_ind])
    )
    values(vmr.gr) <- DataFrame(index[,-(1:2)])
    return(vmr.gr)
  } else {
    start_inds <- sapply(index, function(ix) ix[1])
    end_inds <- sapply(index, function(ix) ix[length(ix)])
    num_cpg <- lengths(index)
    cr.gr <- GRanges(
      seqnames = seqnames(gr)[start_inds],
      ranges = IRanges(start = start(gr)[start_inds],
                       end = end(gr)[end_inds])
    )
    values(cr.gr) <- DataFrame(cr_index = 1:length(index),
                               num_cpg = num_cpg)
    return(cr.gr)
  }
}
