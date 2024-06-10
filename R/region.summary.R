#' @title Compute regional methylation information for individual cells.
#' 
#' @description This function summarize the methylated CpG count and total CpG count
#' per region per cell. 
#'
#' @param region_ranges \code{GRanges} object that contains genomic coordinates
#' of regions of interest.
#' @param SE \code{SummarizedExperiment} object with one (and only one) assay that
#' contains *binary* methylation status of CpG sites in individual cells. In usual
#' analysis workflow (of vmrseq), \code{SE} should be the output of \code{vmrseq::data.pool}.
#' @param sparseNAdrop logical value that represents whether the NA values are
#' droppped in the input \code{SE} object. \code{SE} objects output by
#' \code{vmrseq::data.pool} are NA dropped. See \code{?vmrseq::data.pool}
#' for details about NA-dropped representation.
#' @param verbose logical value that indicates whether progress messages
#' should be printed to stdout. Defaults value is TRUE.
#' @param BPPARAM a \code{BiocParallelParam} object to specify the parallel
#' backend. The default option is \code{BiocParallel::bpparam()} which will
#' automatically creates a cluster appropriate for the operating system.
#' @importFrom GenomicRanges findOverlaps granges
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom S4Vectors subjectHits queryHits
#' @return Returns a \code{SummarizedExperiment} object of dimension 
#' (# regions x # cells). In total the object contains thress assays. Specifically, 
#' `M` and `Cov` represent the number of methylated CpGs and the number of covered 
#' CpGs in each region per cell; and `MF` (stands for methylation fraction) represents 
#' the regional average methylation level computed by `M/Cov`.
#' @export
#'
#' @examples
#' 
#' # load example data
#' toy.se <- HDF5Array::loadHDF5SummarizedExperiment(system.file("extdata", "toy", package = "vmrseq"))
#' data(toy.results)
#' 
#' # run vmrseq.fit
#' regions.se <- region.summary(SE = toy.se, region_ranges = toy.results$vmr.ranges[1:3])
#' regions.se
#' 
region.summary <- function(
    SE,
    region_ranges,
    sparseNAdrop = is_sparse(assays(SE)[[1]]),
    verbose = TRUE,
    BPPARAM = BiocParallel::bpparam()
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
  
  ## Fing the CpG sites that locate within the targeted regions
  hits <- GenomicRanges::findOverlaps(GenomicRanges::granges(SE), region_ranges)

  if (length(hits) > 0) {
    
    idx <- unique(S4Vectors::subjectHits(hits))

    M_mat <- assays(SE)$M_mat

    computeRegionStats <- function(i) { # i th feature/window
      
      mat <- M_mat[S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits)==i], ] %>% 
        matrix(ncol = ncol(SE))
      
      if (!sparseNAdrop) {
        return(list(
          "M" = colSums(mat, na.rm = TRUE),
          "Cov" = colSums(mat >= 0, na.rm = TRUE)
        ))
      } else {
        return(list(
          "M" = round(colSums(mat)),
          "Cov" = colSums(mat > 0)
        ))
      }
    }
    lists <- bplapply(idx, computeRegionStats)
    M <- do.call(
      rbind,
      lapply(lists, function(l) l$M)
    )
    Cov <- do.call(
      rbind,
      lapply(lists, function(l) l$Cov)
    )

    regions.se <- SummarizedExperiment::SummarizedExperiment(
      assays = list("M" = M, "Cov" = Cov, "MF" = M/Cov),
      rowRanges = region_ranges[idx]
    )
  }

  return(regions.se)
}
