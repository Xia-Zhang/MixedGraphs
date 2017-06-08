#' Test function for ADMM Solver
#'
#' @param X the input matrix
#' @param y the responce vector
#' @param method the solver method
#' @param lambda a constant scalar parameter to control the influence of L1-Norm
#' @param ... other parameters
#'
#' @examples
#' X <- matrix(rnorm(500), ncol = 10)
#' y <- rpois(50, 3)
#' glm(y ~ X, family = "poisson")
#' test(cbind(1, X), y, "poisson", 0.0)
#'
test <- function(X, y, method = "Gaussian", lambda = 0.5) {
    .Call('MixedGraphs_test', PACKAGE = 'MixedGraphs', X, y, method, lambda)
}

#' Test function for Newton Solver
#'
#' @param X the input matrix
#' @param y the responce vector
#' @param method the solver method
#' @param lambda a constant scalar parameter to control the influence of L2-Norm
#' @param ... other parameters
#'
#' @examples
#' X <- matrix(rnorm(500), ncol = 10)
#' y <- rpois(50, 3)
#' glm(y ~ X, family = "poisson")
#' testNewton(cbind(1, X), y, "poisson", 0.0)

testNewton <- function(X, y, method = "Gaussian", lambda = 0.5) {
    .Call('MixedGraphs_testNewton', PACKAGE = 'MixedGraphs', X, y, method, lambda)
}