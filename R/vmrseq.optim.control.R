#' @title Auxiliary function as user interface for vmrseq optimization.
#'
#' @description Typically only used when calling vmrseq function with the
#' option \code{control}.
#'
#' @param inits vector of numeric values between 0 and 1 representing initial
#' values of pi_1 shall be taken in optimization algorithm.
#' @param epsilon numeric value representing the convergence upper bound for
#' the algorithm.
#' @param backtrack logical value indicating whether to use backtracking line
#' search to automatically adjust learning rate. Default is TRUE.
#' @param eta a numeric value representing the learning rate in optimization.
#' Default is \code{ifelse(backtrack, 0.05, 0.005)}.
#' @param maxIter positive integer value representing the maximum number of
#' iterations in optimization algorithm.
#'
#' @return the list of arguments for optimization control
#'
#' @export
#'
vmrseq.optim.control <- function(
    inits = c(.2, .5, .8),
    epsilon = 1e-3,
    backtrack = T,
    eta = ifelse(backtrack, 0.05, 0.005),
    maxIter = 100
) {
  return(list(inits = inits, epsilon = epsilon,
              backtrack = backtrack,
              eta = eta, maxIter = maxIter))
}
