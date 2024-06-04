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
HDF5NAdrop2matrix <- function(hdf5_assay) {
  mat <- hdf5_assay |> as("matrix") |> as("sparseMatrix") |> recommenderlab::dropNA2matrix()
  return(mat)
}
