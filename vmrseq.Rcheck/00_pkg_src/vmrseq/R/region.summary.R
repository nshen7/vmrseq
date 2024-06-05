#' @title Compute regional average methylation for individual cells.
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
#' @return Returns a \code{SummarizedExperiment} object that contains the regional
#' average methylation per cell.
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

    SE_M <- assays(SE)$M_mat
    if (sparseNAdrop) SE_M <- SE_M %>% matrix(ncol = ncol(SE)) %>% as("sparseMatrix")

    computeRegionStats <- function(i, type) { # i th feature/window
      mat <- SE_M[S4Vectors::queryHits(hits)[S4Vectors::subjectHits(hits)==i], ]

      if (!sparseNAdrop) {
        if (type == "M") return(colSums(mat, na.rm = T))
        else if (type == "Cov") return(colSums(mat >= 0, na.rm = T))
        else stop("Wrong 'type' value. Either 'Cov' or 'M'.")
      } else {
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
      assays = list("M" = M, "Cov" = Cov),
      rowRanges = region_ranges[idx]
    )
  }

  return(regions.se)
}
