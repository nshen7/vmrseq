#' Title
#'
#' @param gr
#' @param alpha
#' @param maxGap
#' @param minNumCR
#' @param minNumVMR
#' @param minNumLong
#' @param gradient
#' @param tp
#' @param control
#' @param verbose
#' @param BPPARAM
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
    maxNumMerge = 0,
    minNumLong = 0,
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
  cutoff <- computeVarCutoff(
    alpha = alpha,
    meth = gr$meth,
    total = gr$total
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
    pct_incr <- round(sum(lengths(CRI))/length(SE)*100, 2)
    message("...Finished calling candidate regions - found ", length(CRI),
            " candidate regions in total.
  ...", pct_incr,
            "% QC-passed sites are called to be in candidate regions.")

    # Add summary stats (smoothed var and CR index) into output
    cr_index <- rep(NA, length(SE))
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
  ...", round(sum(vmr.df$end_ind-vmr.df$start_ind+1) / length(SE) * 100, 2),
              "% QC-passed sites are called to be in VMRs.")
    }

    # Formatting function output
    vmr_index <- rep(NA, length(SE))
    vmr_index[unlist(VMRI)] <- rep.int(1:length(VMRI), lengths(VMRI))
    values(gr) <- cbind(values(gr), vmr_index)

    vmr.gr <- indexToGranges(gr = gr, Indexes = VMRI)
    values(vmr.gr) <- cbind(values(vmr.gr), vmr.df[, -c(3,5)])

    cr.gr <- indexToGranges(gr = gr, Indexes = CRI)

    return(list(gr = gr, vmr.ranges = vmr.gr, cr.ranges = cr.gr, alpha = alpha, var_cutoff = cutoff))
  }
}
