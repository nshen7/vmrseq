#' Title
#'
#' @param gr
#' @param cutoff
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
#' @import dplyr
#' @import GenomicRanges
#'
#' @return
#' @export
#'
#' @examples
vmrseq.fit <- function(
    gr,
    cutoff = 0.1,
    maxGap = 1000,
    minNumCR = 5, minNumVMR = 5,
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

  cr_index <- rep(NA, length(SE))
  cr_index[unlist(CRI)] <- rep.int(1:length(CRI), lengths(CRI))

  # Add summary stats (smoothed var and CR index) into output
  values(gr)$cr_index <- cr_index

  # Percentage of sites in CRs
  pct_incr <- round(sum(lengths(CRI))/length(SE)*100, 2)
  if (is.null(CRI)) {
    message("...No candidate regions pass the cutoff")
    return(NULL)
  } else {
    message("...Finished calling candidate regions - found ", length(CRI),
            " candidate regions in total.
  ...", pct_incr,
            "% QC-passed sites are called to be in candidate regions.")
  }


  message(
    "Step 2: Detecting VMRs",
    # ifelse(penalty > 0, yes = " with", no = " without"),
    # " penalty",
    # ifelse(penalty > 0, yes = paste0("=", round(penalty, 2)), no = ""),
    "..."
  )
  t1 <- proc.time()

  vmrs.df <- searchVMR( # data frame of VMR information
    gr = gr,
    CRI = CRI,
    # penalty = penalty,
    maxGap = maxGap,
    minNumVMR = minNumVMR,
    minNumLong = minNumLong,
    maxNumMerge = maxNumMerge,
    tp = tp,
    gradient = gradient,
    control = control,
    verbose = verbose,
    parallel = parallel
  )
  VMRI <- lapply( # list of indices of VMRs
    1:nrow(vmrs.df),
    function(i) vmrs.df$start_ind[i]:vmrs.df$end_ind[i]
  )

  t2 <- proc.time()

  if (nrow(vmrs.df) == 0) {
    message("No VMR detected.")
    return(NULL)
  } else {
    message("...Finished detecting VMRs - took ",
            round((t2 - t1)[3]/60, 2), " min and ",
            nrow(vmrs.df), " VMRs found in total.
  ...", round(sum(vmrs.df$end_ind-vmrs.df$start_ind+1) / length(SE) * 100, 2),
            "% QC-passed sites are called to be in VMRs.")
  }

  vmr.gr <- indexToGranges(gr = gr, Indexes = VMRI)
  cr.gr <- indexToGranges(gr = gr, Indexes = CRI)

  vmr_index <- rep(NA, length(SE))
  vmr_index[unlist(VMRI)] <- rep.int(1:length(VMRI), lengths(VMRI))
  values(gr)$vmr_index <- vmr_index

  return(list(gr = gr, vmr.ranges = vmr.gr, cr.ranges = cr.gr))
}
