#' @title Compute regional methylation information for individual cells.
#' 
#' @description This function summarize the methylated CpG count and total CpG count
#' per region per cell.#' 
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
#'
region.summary <- function(
    SE,
    region_ranges,
    sparseNAdrop = is_sparse(assays(SE)[[1]])
) {

  hits <- GenomicRanges::findOverlaps(GenomicRanges::granges(SE), region_ranges)

  if (length(hits) > 0) {
    
    idx <- unique(S4Vectors::subjectHits(hits))

    M_mat <- assays(SE)$M_mat

    computeRegionStats <- function(i, type) { # i th feature/window
      
      mat <- M_mat[S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits)==i], ] %>% 
        matrix(ncol = ncol(SE))
      
      if (!sparseNAdrop) {
        # mat <- M_mat[S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits)==i], ]
        if (type == "M") return(colSums(mat, na.rm = T))
        else if (type == "Cov") return(colSums(mat >= 0, na.rm = T))
        else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
      } else {
        # mat <- M_mat[S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits)==i], ] %>% 
        #   matrix(ncol = ncol(SE)) %>% 
        #   as("sparseMatrix")
        if (type == "M") return(round(colSums(mat)))
        else if (type == "Cov") return(colSums(mat > 0))
        else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
      }
    }
    M <- do.call(
      rbind,
      bplapply(idx, computeRegionStats, type = "M")
    )
    Cov <- do.call(
      rbind,
      bplapply(idx, computeRegionStats, type = "Cov")
    )

    regions.se <- SummarizedExperiment::SummarizedExperiment(
      assays = list("M" = M, "Cov" = Cov, "MF" = M/Cov),
      rowRanges = region_ranges[idx]
    )
  }

  return(regions.se)
}
