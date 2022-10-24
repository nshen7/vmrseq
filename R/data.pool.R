#' Pool single-cell file together into an HDF5-based SummarizedExperiment
#' object (with or w/o sparse matrix representation)
#'
#' @param cellFiles Vector of character strings indicating single-cell file
#' paths you wish to pool into SE object(s). Cell files should be in BED-like
#' format, where the first 5 columns in each file must be: <chr>, <pos>, <strand>,
#' <meth_read>, <total_read>, in strict order. Each row contains info of one CpG.
#' Strand is assumed to be properly filpped. Cell files can be zipped ones in .gz
#' format.
#' @param writeDir A single character string indicating a folder directory where
#' you wish to store the processed SE object(s).
#' @param selectChrs Either NULL or a vector of character strings. If
#' \code{selectChrs = NULL}, all chromosomes in \code{cellFiles} will be processed;
#' otherwise, only chromosomes listed in \code{selectChrs} will be processed.
#' @param sepChrs Logical value indicating whether or not to separate chromosomes
#' in different SE objects during processing and saving to disk. Default is TRUE.
#' @param sparseNAdrop Logical value indicating whether or not to use NA-dropped
#' sparseMatrix representation. 'NA-dropped' means replacing NAs as zeros and
#' then represent 0 as a very small positive value close to 0 so that the data
#' matrix is converted to a traditional sparseMatrix form with 0 as default
#' entry value. (see \code{sparseNAMatrix-class} from package \code{recommenderlab}).
#' Using NA-dropped sparseMatrix representation helps to save RAM during data
#' processing, as well as storage space saving to disk. Default value is TRUE.
#'
#' @return
#' @export
#'
#' @examples
#'

data.pool <- function(cellFiles,
                      writeDir,
                      selectChrs = NULL,
                      sepChrs = TRUE,
                      sparseNAdrop = TRUE) {

  # TODO: making checks on input data format

  if (!sepChrs) stop("Sorry, `sepChrs=FALSE` is not inplemented yet!")

}
