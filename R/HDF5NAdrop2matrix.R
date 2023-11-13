#' Coerce a NA-dropped HDF5Matrix back to matrix form
#'
#' @param HDF5Matrix an HDF5Matrix (e.g. assay output from vmrseq::data.pool)
#' with dropped NA values
#'
#' @importFrom recommenderlab dropNA2matrix
#' @return
#' @export
#'
#' @examples
HDF5NAdrop2matrix <- function(hdf5_matrix) {
  mat <- hdf5_matrix %>% as("sparseMatrix") %>% recommenderlab::dropNA2matrix()
  return(mat)
}