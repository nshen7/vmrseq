% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vmrseq.smooth.R
\name{vmrseq.smooth}
\alias{vmrseq.smooth}
\title{Smoothing on single-cell bisulfite sequencing data for the purpose of
constructing candidate regions.}
\usage{
vmrseq.smooth(
  SE,
  bpWindow = 2000,
  sparseNAdrop = is_sparse(assays(SE)[[1]]),
  verbose = TRUE,
  BPPARAM = bpparam()
)
}
\arguments{
\item{SE}{\code{SummarizedExperiment} object with one (and only one) assay that
contains *binary* methylation status of CpG sites in individual cells. We
recommend using output by \code{vmrseq::data.pool} (i.e., an NA-dropped
HDF5-based SummarizedExperiment object) to prevent running out of memory.}

\item{bpWindow}{positive integer that represents the width (in bp) of
smoothing window. Default value is 2000.}

\item{sparseNAdrop}{logical value that represents whether the NA values are
droppped in the input \code{SE} object. \code{SE} objects output by
\code{vmrseq::data.pool} are NA dropped. See \code{?vmrseq::data.pool}
for details about NA-dropped representation.}

\item{verbose}{logical value that indicates whether progress messages
should be printed to stdout. Defaults value is TRUE.}

\item{BPPARAM}{a \code{BiocParallelParam} object to specify the parallel
backend. The default option is \code{BiocParallel::bpparam()} which will
automatically creates a cluster appropriate for the operating system.}
}
\value{
a \code{GRanges} object that contains the result of smoothing.
 The object retains genomic coordinates (chr, start, end) of input CpG
 sites, in the same order as in the input \code{SE} object. Three
 column are added (on top of original metadata columns for the CpG sites in
 \code{SE}, if any):
 1. meth: methylated cell count of the CpG
 2. total: total (non-missing) cell count of the CpG
 3. var: variance computed based on individual-cell smoothed relative methylation levels.
}
\description{
\code{vmrseq.smooth} takes a \code{SummarizedExperiment} object
 with information of methylation level of individual cells as input, and
 perform a kernel smoother to ‘relative’ methylation levels of individual
 cells prior to constructing candidate regions. Purpose of the smoothing is
 to adjust for uneven coverage biases and borrow information from nearby sites.
 See manuscript for detailed description.
}
\examples{
# load example data
toy.se <- HDF5Array::loadHDF5SummarizedExperiment(system.file("extdata", "toy", package = "vmrseq"))

# preprocessing
total <- DelayedArray::rowSums(SummarizedExperiment::assays(toy.se)$M_mat > 0)
toy.se <- subset(toy.se, total >= 3)

# run vmrseq.smooth
toy.gr <- vmrseq.smooth(toy.se)
toy.gr

}
\seealso{
\code{\link{data.pool}}, \code{\link{vmrseq.fit}}
}
