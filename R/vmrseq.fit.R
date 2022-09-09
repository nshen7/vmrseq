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

}
