#' Coerce a NA-dropped HDF5Matrix back to matrix form
#'
#' @param HDF5Matrix an HDF5Matrix (e.g. an assay from vmrseq::data.pool output)
#' with dropped NA values
#'
#' @return
#' @export
#'
#' @examples
HDF5NAdrop2matrix <- function(hdf5_matrix) {
  mat <- hdf5_matrix %>% as("sparseMatrix") %>% recommenderlab::dropNA2matrix()
  return(mat)
}
