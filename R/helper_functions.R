callCandidRegion <- function(SE,
                             qVar,
                             maxGap, minNumCR,
                             bpWindow,
                             acrsSmooth, bpSpan, minInSpan,
                             maxNumMerge,
                             verbose,
                             parallel) {

  # gr <- granges(SE)
  # M <- assays(SE)[[1]]
  #
  # # Compute methylated fraction of cells for individual sites
  # gr$mean_meth <- gr$meth / gr$total
  #
  # # Apply smoother and compute variance on each chromosome serially
  # chrs <- as.character(unique(seqnames(gr)))
  # fit <- NULL; var <- NULL
  # for (chromosome in chrs) {
  #
  #   if (verbose) message(
  #     "...Chromosome ",
  #     paste(chromosome, collapse = ", "), ": ",
  #     appendLF = FALSE
  #   )
  #
  #   t1 <- proc.time()
  #
  #   # Subset one chromosome from gr and M
  #   gr_chr <- subset(gr, seqnames(gr) == chromosome)
  #   if (length(gr_chr) < minNumCR){
  #     message("No candidates found.")
  #     next
  #   }
  #   M_chr <- M[seqnames(gr) == chromosome, ]
  #
  #   ## Locfit smooth on fractional methylation
  #   if (acrsSmooth) {
  #     weights <- gr_chr$total
  #     fit_chr <- smoothMF(x = start(gr_chr), y = gr_chr$mf,
  #                         chr = chromosome,
  #                         weights = weights,
  #                         maxGap = maxGap, minNumCR = minNumCR,
  #                         minInSpan = minInSpan, bpSpan = bpSpan,
  #                         verbose = verbose,
  #                         parallel = parallel)
  #     fit_chr$fitted <- fit_chr$fitted %>% pmax(0) %>% pmin(1)
  #   } else {
  #     fit_chr <- data.frame(fitted = rep(NA, length(gr_chr)), is_smooth = FALSE)
  #   }
  #   # Keep the raw mf if not smoothed
  #   ind <- which(!fit_chr$is_smooth)
  #   if (length(ind) > 0) fit_chr$fitted[ind] <- gr_chr$mf[ind]
  #   # Concatenate fit_chr from all chromosomes
  #   fit <- rbind(fit, fit_chr)
  #   if (verbose) message("Smoothed, ", appendLF = FALSE)
  #
  #
  #   ## Compute variance relative to smoothed MF
  #   var_chr <- computeVar(gr_chr, M_chr, fit_chr,
  #                               bpWindow,
  #                               parallel)
  #   # Concatenate fit_chr from all chromosomes
  #   var <- c(var, var_chr)
  #
  #   if (verbose) {
  #     t2 <- proc.time()
  #     message("variance computed (", round((t2 - t1)[3]/60, 2), " min). ")
  #   }
  # }
  # colnames(fit) <- c("smoothed_mf", "is_smooth")

  # Compute cutoff based on qVar
  cutoff <- quantile(var, prob = 1-qVar)
  mes <- "...Calling candidate regions with cutoff of %.2f on variance."
  message(sprintf(mes, cutoff))

  # Call candidate regions
  cluster <- bumphunter::clusterMaker(chr = seqnames(gr),
                                      pos = start(gr),
                                      maxGap = maxGap,
                                      assumeSorted = TRUE)
  Indexes <- bumphunter::getSegments(x = var, f = cluster,
                                     cutoff = cutoff,
                                     assumeSorted = TRUE,
                                     verbose = FALSE)

  upIndex <- Indexes$upIndex
  if (length(upIndex) == 0) upIndex <- NULL

  if (!is.null(upIndex)) {
    # Only keep candidate regions with more than `minNumCR` CpGs
    CRI <- upIndex[lengths(upIndex) >= minNumCR]
    return(list(CRI = CRI, mean_meth_smooth = fit, var = var))
  } else {
    return(list(CRI = NULL, mean_meth_smooth = fit, var = var))
  }

} # end of function `callCandidRegion`



#' Title
#'
#' @param x
#' @param y
#' @param weights
#' @param chr
#' @param minInSpan
#' @param bpSpan
#' @param verbose
#' @param parallel
#'
#' @return
#' @export
#'
#' @examples
smoothMF <- function(x, y, weights, chr,
                     minInSpan, bpSpan,
                     verbose,
                     parallel) {


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
      stop("y (mean_meth) is missing")
    if (is.null(xi))
      stop("x (pos) is missing")
    if (is.null(clusteri))
      stop("cluster is missing")

    if (length(idx) >= minInSpan) {
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
    return(data.frame(fitted = as.vector(yi), is_smooth = is_smooth))
  } # end of function `locfitByCluster`

  chr_len <- rep(chr, each = length(x))
  clusterC <- bumphunter::clusterMaker(factor(chr_len, levels=unique(chr_len)),
                                       x,
                                       assumeSorted = TRUE,
                                       maxGap = 2*bpSpan)

  Indexes <- split(seq(along = clusterC), clusterC)

  if (parallel) {
    ret <- do.call(rbind, bplapply(Indexes, function(idx) locfitByCluster(idx)))
  } else {
    ret <- do.call(rbind, lapply(Indexes, function(idx) locfitByCluster(idx)))
  }

  return(ret) # data.frame with columns: 'fitted' (numeric), 'is_smooth' (logical)

} # end of function `smoothMF`


#' Title
#'
#' @param gr
#' @param M numeric matrix of binary single-cell methylation status. Row number
#' should be equal to the number of CpG sites and column number should be equal
#' to the number of cells.
#' @param mean_meth numeric vector of across-cell mean methylation. Length should
#' be equal to the number of CpG sites.
#' @param bpWindow
#' @param parallel
#'
#' @return
#' @export
#'
#' @examples
computeVar <- function(gr,
                       M,
                       mean_meth,
                       bpWindow,
                       parallel) {

  varByCluster <- function(idx) {

    gr_i <- gr[idx, ]
    mean_meth_i <- mean_meth[idx, ]
    if (length(idx) > 1) {
      M_i <- M[idx, ]
    } else {
      M_i <- t(M[idx, ])
    }

    # Obtain window ranges of each site
    wds_i <- GRanges(
      seqnames = seqnames(gr_i),
      ranges = IRanges(start = start(gr_i) - round(bpWindow/2),
                       end = end(gr_i) + round(bpWindow/2))
    )
    # Find indices in gr for each window
    hits <- findOverlaps(wds_i, gr_i)
    wds_inds <- split(subjectHits(hits), queryHits(hits))

    # Matrix of 'relative methylation' or 'methylation residuals'
    rel_M_i <- as.matrix(M_i - mean_meth_i)

    varByWindow <- function(j) {
      ind <- wds_inds[[j]]
      if (length(ind) > 1) {
        # # local constant smoothing with Gaussian kernel
        # rel_M_i_j <- rel_M_i[ind, ]
        # gr_i_j <- gr_i[ind]
        # d <- (start(gr_i_j) - start(gr_i[j])) / (bpWindow/2)
        # # w <- (1 - abs(d)^3)^3 # tricube kernel
        # w <- exp(-(2.5*d)^2/2) # gaussian kernel
        # w_mat <- matrix(w, nrow = length(d), ncol = ncol(rel_M_i), byrow = FALSE)
        # w_mat[which(is.na(rel_M_i_j))] <- NA
        # wd_mean_meths <- colSums(rel_M_i_j*w_mat, na.rm = TRUE)/colSums(w_mat, na.rm = TRUE)

        # local constant smoothing with box kernel
        wd_mean_meths <- colMeans(rel_M_i[ind, ], na.rm = TRUE)
      } else {
        wd_mean_meths <- rel_M_i[ind, ]
      }
      wd_var <- var(wd_mean_meths, na.rm = TRUE)
      return(wd_var)
    }

    var_i <- do.call(c, lapply(1:length(wds_inds), function(j) varByWindow(j)))
    # mean_meth_sm_i <- rowMeans(rel_M_i, na.rm = TRUE)

    return(var_i)
  } # end of function 'varByCluster'

  clusterC <- bumphunter::clusterMaker(seqnames(gr),
                                       start(gr),
                                       assumeSorted = TRUE,
                                       maxGap = bpWindow + 1)
  Indexes <- split(seq(along = clusterC), clusterC)

  if (parallel) {
    var <- do.call(c, bplapply(Indexes, varByCluster)) %>% unname
  } else {
    var <- do.call(c, lapply(Indexes, varByCluster)) %>% unname
  }

  return(var)
} # end of function 'computeVar'


searchVMR <- function(gr,
                      CRI,
                      penalty,
                      maxGap, minNumVMR,
                      tp,
                      maxNumMerge,
                      minNumLong,
                      gradient,
                      control,
                      verbose,
                      parallel) {

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

    res_1g <- .solve1Grp(pos = pos, totals = totals, meths = meths,
                         tp = tp,
                         METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)

    if (length(ix) >= minNumLong) { # long region
      res_2g <- .solve2Grp(gradient = gradient,
                           pos = pos, totals = totals, meths = meths,
                           tp = tp,
                           inits = control$inits, epsilon = control$epsilon,
                           backtrack = control$backtrack,
                           eta = control$eta, max_iter = control$maxIter,
                           CHOICEARRAY = CHOICEARRAY,
                           METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
    } else { # short region
      res_2g <- .solve2Grp(gradient = gradient,
                           pos = pos, totals = totals, meths = meths,
                           tp = tp,
                           inits = mean(meths/totals), epsilon = control$epsilon,
                           backtrack = control$backtrack,
                           eta = control$eta, max_iter = control$maxIter,
                           CHOICEARRAY = CHOICEARRAY,
                           METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)
    }

    # if (res_2g$loglik > res_1g$loglik + penalty) {
    if (TRUE) {
      vmr_inds <- .callVMR(
        state_seq_2g = res_2g$vit_path[, 1:2],
        min_n = minNumVMR,
        max_n_merge = maxNumMerge
      )
      if (is.null(vmr_inds)) return(NULL)
      else return(data.frame(vmr_inds + ix[1] - 1,
                             # cr_index = i,
                             vmr_num_cpg = vmr_inds$end_ind - vmr_inds$start_ind + 1,
                             pi = res_2g$optim_pi_1,
                             n_iter = res_2g$n_iter,
                             loglik_diff = res_2g$loglik - res_1g$loglik))
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
