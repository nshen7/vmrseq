#' Coerce a NA-dropped HDF5 matrix back to regular matrix form where NAs are not
#' dropped to 0.
#'
#' @param hdf5_assay an DelayedMatrix object with dropped NA values (e.g. an
#' assay from HDF5SummarizedExpriment object saved by vmrseq::data.pool)
#'
#' @importFrom recommenderlab dropNA2matrix
#' @return Returns a matrix.

#' @export
#' 
#' @examples
#' # load example data
#' toy.se <- HDF5Array::loadHDF5SummarizedExperiment(system.file("extdata", "toy", package = "vmrseq"))
#' 
#' # run the function
#' HDF5NAdrop2matrix(SummarizedExperiment::assays(toy.se)$M_mat)
#' 
#'
HDF5NAdrop2matrix <- function(hdf5_assay) {
  mat <- hdf5_assay |> as("matrix") |> as("sparseMatrix") |> recommenderlab::dropNA2matrix()
  return(mat)
}
