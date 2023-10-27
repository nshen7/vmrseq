#' Construct candidate regions (CRs) by taking groups of consecutive loci that exceed
#' threshold on the variance of smoothed relative methylation levels and detect
#' variably methylated regions (VMRs) by optimizing a hidden Markov model (HMM).
#'
#' @param gr \code{GRanges} object output by \code{vmrseq::vmrseq.smooth},
#' containing genomic coordinates (chr, start, end) and summarized information
#' (meth, total, var) of CpG sites in the input dataset.
#' @param alpha positive scalar value between 0 and 1 that represents the
#' designated significance level for determining variance threshold of candidate
#' regions construction. The variance threshold is determined by taking 1-alpha
#' quantile value of an approximate null distribution of variance simulated
#' from the beta priors of emission probability in the hidden Markov model.
#' Default value of alpha is 0.05.
#' @param maxGap integer value representing maximum number of base pairs in
#' between neighboring CpGs to be included in the same VMR. Default value is
#' 2000 bp.
#' @param minNumCR positive integer value representing the minimum number of
#' CpG sites within a candidate region. Default value is 5.
#' @param minNumVMR positive integer value representing the minimum number of
#' CpG sites within a variably methylated region. Default value is 5.
#' @param gradient logical value indicating whether exponentiated gradient
#' descent shall be applied to update prevalence parameter. Default is TRUE. If
#' set as FALSE, initial values (i.e., value of \code{inits} arguments in
#' \code{vmrsqe::vmrseq.optim.control}, can be set up in the \code{control}
#' argument of this function) are used as prevalence parameter for decoding
#' hidden states.
#' @param tp a `transitProbs-class` object that contains the transition
#' probability distribution used for HMM optimization. Default value is
#' transition probability \code{vmrseq:::tp0} built in the package that was
#' previously trained on mouse brain cells. See manuscript for training
#' procedure and data source.
#' @param control list of miscellaneous parameters used to control optimization
#' of the HMM model. Default is output of \code{vmrseq::vmrseq.optim.control()}.
#' Can be changed by tweaking arguments in function \code{vmrseq::vmrseq.optim.control()}.
#' @param verbose logical value that indicates whether progress messages
#' should be printed to stdout. Defaults value is TRUE.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#' @importFrom BiocParallel bplapply register MulticoreParam bpparam
#' @importFrom bumphunter clusterMaker getSegments
#' @importFrom gamlss.dist rBEZI rBE
#' @import dplyr
#' @import GenomicRanges
#'
#' @return
#' @export
#'
#' @examples
vmrseq.fit <- function(
    gr,
    alpha = 0.05,
    maxGap = 2000,
    minNumCR = 5, minNumVMR = 5,
    gradient = TRUE,
    tp = NULL,
    control = vmrseq.optim.control(),
    verbose = TRUE,
    BPPARAM = bpparam()
) {

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

  # Compute cutoff from beta priors
  pars <- getPriorParams(gr$total)
  cutoff <- computeVarCutoff(
    alpha = alpha,
    meth = gr$meth,
    total = gr$total,
    pars_u =  pars$pars_u,
    pars_m = pars$pars_m
  )

  # Bumphunt candidate regions
  message("Step 1: Detecting candidate regions...")
  CRI <- callCandidRegion(
    gr = gr,
    cutoff = cutoff,
    maxGap = maxGap,
    minNumCR = minNumCR,
    bpWindow = bpWindow,
    verbose = verbose,
    parallel = parallel
  ) # Outputs list of index vectors. Each list element contains indices in of a CR `gr`.

  if (is.null(CRI)) { # End the function if no CR detected
    message("...No candidate regions pass the cutoff")
    return(NULL)
  } else {
    pct_incr <- round(sum(lengths(CRI))/length(gr)*100, 2)
    message("...Finished calling candidate regions - found ", length(CRI),
            " candidate regions in total.
  ...", pct_incr,
            "% QC-passed sites are called to be in candidate regions.")

    # Add summary stats (smoothed var and CR index) into output
    cr_index <- rep(NA, length(gr))
    cr_index[unlist(CRI)] <- rep.int(1:length(CRI), lengths(CRI))
    values(gr)$cr_index <- cr_index

    # Starting detecting VMRs
    message(
      "Step 2: Detecting VMRs..."
    )
    t1 <- proc.time()

    vmr.df <- searchVMR( # data frame of VMR information
      gr = gr,
      CRI = CRI,
      # penalty = penalty,
      # maxGap = maxGap,
      minNumVMR = minNumVMR,
      # minNumLong = minNumLong,
      # maxNumMerge = maxNumMerge,
      tp = tp,
      gradient = gradient,
      control = control,
      verbose = verbose,
      parallel = parallel
    )
    VMRI <- lapply( # list of indices of VMRs
      1:nrow(vmr.df),
      function(i) vmr.df$start_ind[i]:vmr.df$end_ind[i]
    )

    t2 <- proc.time()

    if (nrow(vmr.df) == 0) {
      message("No VMR detected.")
      return(NULL)
    } else {
      message("...Finished detecting VMRs - took ",
              round((t2 - t1)[3]/60, 2), " min and ",
              nrow(vmr.df), " VMRs found in total.
  ...", round(sum(vmr.df$end_ind-vmr.df$start_ind+1) / length(gr) * 100, 2),
              "% QC-passed sites are called to be in VMRs.")
    }

    # Formatting function output
    vmr_index <- rep(NA, length(gr))
    vmr_index[unlist(VMRI)] <- rep.int(1:length(VMRI), lengths(VMRI))
    values(gr) <- cbind(values(gr), vmr_index)

    vmr.gr <- indexToGranges(gr = gr, Indexes = VMRI)
    values(vmr.gr) <- cbind(values(vmr.gr), vmr.df[, -c(3,5)])

    cr.gr <- indexToGranges(gr = gr, Indexes = CRI)

    return(list(gr = gr, vmr.ranges = vmr.gr, cr.ranges = cr.gr, alpha = alpha, var_cutoff = cutoff, bb_params = pars))
  }
}
