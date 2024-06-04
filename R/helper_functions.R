#' cell_1
#'
#' This dataset is an example of a single-cell file that can be input to
#' \code{vmrseq::data.pool} function. It contains the first 100 CpG sites in
#' a cell from mouse frontal cortext dataset published by Luo et al. (2017).
#'
#' @docType data
#' @name cell_1
#' @usage data(cell_1)
#' @format A data frame with 100 rows and 5 variables (no column names):
#' \describe{
#'   \item{V1}{Chromosome}
#'   \item{V2}{Genomic coordinate}
#'   \item{V3}{Strand information}
#'   \item{V4}{Number of methylated reads}
#'   \item{V5}{Number of reads in total}
#' }
#' @references Luo, Chongyuan et al. \emph{Single-cell methylomes identify neuronal
#' subtypes and regulatory elements in mammalian cortex.}. Science (New York, N.Y.)
#' vol. 357,6351 (2017): 600-604.
#' @examples
#' data(cell_1)
#' cell_1
#'
"cell_1"



#' Kernel smoothing function
#'
#' @param x vector of the x values
#' @param y vector of the y values
#' @param weights vector of the weight associated with each data point
#' @param chr string of the chromosome name
#' @param minInSpan minimum number of sites in the span
#' @param bpSpan base pair of the span
#' @param verbose logical value that indicates whether progress messages
#' should be printed to stdout
#' @param parallel logical value that indicates whether function should be
#' run in parallel

smoothMF <- function(x, y,
                     weights, chr,
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
      stop("y (meanMeth) is missing")
    if (is.null(xi))
      stop("x (pos) is missing")
    if (is.null(clusteri))
      stop("cluster is missing")

    df <- data.frame(posi = xi, yi = yi, weightsi = weightsi)

    # balance minInSpan and bpSpan
    nn <- minInSpan / length(idx)
    fit <- locfit(yi ~ lp(posi, nn = nn, h = bpSpan),
                  data = df, weights = weightsi, family = "gaussian",
                  maxk = 10000)
    yi <- fitted(fit)

    return(as.vector(yi))

  } # end of function `locfitByCluster`

  chr_len <- rep(chr, each = length(x))
  clusterC <- bumphunter::clusterMaker(factor(chr_len, levels=unique(chr_len)),
                                       x,
                                       assumeSorted = TRUE,
                                       maxGap = 2*bpSpan)

  Indexes <- split(seq(along = clusterC), clusterC)

  if (parallel) {
    ret <- do.call(c, bplapply(Indexes, function(idx) locfitByCluster(idx)))
  } else {
    ret <- do.call(c, lapply(Indexes, function(idx) locfitByCluster(idx)))
  }

  return(ret) # data.frame with columns: 'fitted' (numeric), 'is_smooth' (logical)

} # end of function `smoothMF`


#' Compute the smoothed variance
#'
#' @param gr same as in \code{vmrseq.smooth}
#' @param M numeric matrix of binary single-cell methylation status. Row number
#' should be equal to the number of CpG sites and column number should be equal
#' to the number of cells.
#' @param meanMeth numeric vector of across-cell mean methylation. Length should
#' be equal to the number of CpG sites.
#' @param bpWindow same as in \code{vmrseq.smooth}
#' @param sparseNAdrop same as in \code{vmrseq.smooth}
#' @param parallel logical value that indicates whether function should be
#' run in parallel

computeVar <- function(gr, M,
                       meanMeth,
                       bpWindow,
                       sparseNAdrop,
                       parallel) {

  varByCluster <- function(idx) {

    gr_i <- gr[idx, ]
    meanMeth_i <- meanMeth[idx]

    if (!sparseNAdrop) {
      M_i <- matrix(M[idx, ], nrow = length(idx))
    } else {
      M_i <- matrix(M[idx, ] %>% as("sparseMatrix") %>% dropNA2matrix(), nrow = length(idx))
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
    rel_M_i <- M_i - meanMeth_i

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
        # wd_meanMeths <- colSums(rel_M_i_j*w_mat, na.rm = TRUE)/colSums(w_mat, na.rm = TRUE)

        # local constant smoothing with box kernel
        wd_meanMeths <- colMeans(rel_M_i[ind, ], na.rm = TRUE)
      } else {
        wd_meanMeths <- rel_M_i[ind, ]
      }
      wd_var <- var(wd_meanMeths, na.rm = TRUE)
      return(wd_var)
    }

    var_i <- do.call(c, lapply(1:length(wds_inds), function(j) varByWindow(j)))

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

getPriorParams <- function(total) {
  pars_u <- .priorParams(median(total), type = "u")
  pars_m <- .priorParams(median(total), type = "m")
  return(list(pars_u = pars_u, pars_m = pars_m))
}

#' Compute cutoff on smoothed variance for determining candidate regions
#'
#' @param alpha level of significance
#' @param meth vector of meth read counts
#' @param pars_u parameters used in the ZIBB distribution for unmethylated grouping
#' @param pars_m parameters used in the BB distribution for methylated grouping
#' @param n number of simulations
#' @param total vector of total read counts
#'
#' @importFrom gamlss.dist rBEZI rBE

computeVarCutoff <- function(alpha, meth, total, pars_u, pars_m, n = 100000) {
  mf <- meth / total
  prob <- c(sum(mf < 0.4), sum(mf > 0.6)) / sum(mf < 0.4 | mf > 0.6)
  null_p_u <- gamlss.dist::rBEZI(n = round(n*prob[1]), mu = pars_u['mu'], sigma = pars_u['sigma'], nu = pars_u['nu'])
  null_p_m <- gamlss.dist::rBE(n = round(n*prob[2]), mu = pars_m['mu'], sigma = pars_m['sigma'])
  null_p <- c(null_p_u, null_p_m)
  null_var <- null_p * (1-null_p)
  return(quantile(null_var, 1-alpha))
}


callCandidRegion <- function(gr,
                             cutoff,
                             maxGap,
                             minNumCR,
                             bpWindow,
                             verbose,
                             parallel
) {

  ## Compute cutoff based on qVar
  mes <- "...Calling candidate regions with cutoff of %.3f on variance."
  message(sprintf(mes, cutoff))

  ## Remove rows with var==NA
  gr <- subset(gr, !is.na(gr$var))

  ## Call candidate regions
  cluster <- bumphunter::clusterMaker(chr = seqnames(gr),
                                      pos = start(gr),
                                      maxGap = maxGap,
                                      assumeSorted = TRUE)
  Indexes <- bumphunter::getSegments(x = gr$var, f = cluster,
                                     cutoff = cutoff,
                                     assumeSorted = TRUE,
                                     verbose = FALSE)

  upIndex <- Indexes$upIndex
  CRI <- upIndex[lengths(upIndex) >= minNumCR]
  if (length(CRI) == 0) CRI <- NULL

  return(CRI)
} # end of function `callCandidRegion`


searchVMR <- function(gr,
                      CRI,
                      minNumVMR,
                      tp,
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

    res_2g <- .solve2Grp(gradient = gradient,
                         pos = pos, totals = totals, meths = meths,
                         tp = tp,
                         inits = control$inits, epsilon = control$epsilon,
                         backtrack = control$backtrack,
                         eta = control$eta, max_iter = control$maxIter,
                         CHOICEARRAY = CHOICEARRAY,
                         METHARRAY = METHARRAY, UNMETHARRAY = UNMETHARRAY)

    if (res_2g$loglik > res_1g$loglik) {
      vmr_inds <- .callVMR(
        state_seq_2g = res_2g$vit_path[, 1:2],
        min_n = minNumVMR
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
    vmrs.df <- do.call(
      rbind,
      bplapply(1:length(CRI), function(i) searchVMRbyRegion(i))
    )
  } else {
    vmrs.df <- do.call(
      rbind,
      lapply(1:length(CRI), function(i) searchVMRbyRegion(i))
    )
  }

  # return  a data.frame with columns:
  return(vmrs.df)

} # end of function `searchVMR`


indexToGranges <- function(gr, Indexes) {
  start_inds <- sapply(Indexes, function(ix) ix[1])
  end_inds <- sapply(Indexes, function(ix) ix[length(ix)])
  num_cpg <- lengths(Indexes)
  output.gr <- GRanges(
    seqnames = seqnames(gr)[start_inds],
    ranges = IRanges(start = start(gr)[start_inds],
                     end = end(gr)[end_inds])
  )
  values(output.gr) <- DataFrame(num_cpg = num_cpg)
  return(output.gr)
}
