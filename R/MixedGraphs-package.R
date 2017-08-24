#' The package overview
#' 
#' @description This package allows for Gaussian, Logistic and Poisson regression with l1 or l2 penalty. And it can also do fitting and visualizing mixed graphical models.
#'  
#' @name MixedGraphs-package
#' @aliases MixedGraphs-package MixedGraphs
#' @docType package
#' @author Xia Zhang
#' 
#' Maintainer: Xia Zhang <zhangxia9403@@gmail.com>
#' @keywords package
#' @examples
#'
#' ## regressions
#' X <- matrix(rnorm(10 * 50), 10, 50)
#' y <- rnorm(10)
#' glmLasso(X, y, lambda = 0.5, family = "gaussian", support_stability = 10)
#' glmRidge(X, y, lambda = 0.5, family = "gaussian", thresh = 0.005)
#'
#' ## BRAIL
#' X <- lapply(1:2, function(x) {matrix(rnorm(5 * 10), 5, 10)})
#' y <- rnorm(5)
#' BRAIL(X, y, family = "gaussian", tau = 0.8, B = 20, doPar = TRUE)
#'
#' ## MixedGraph fitting and plotting
#' X <- lapply(1 : 2, function(x){matrix(rnorm(12), nrow = 4)})
#' crf_structure = matrix(c(1, 0, 1, 1), 2, 2)
#' brail_control <- list(B = 5, tau = 0.6)
#' G <- MixedGraph(X, crf_structure, brail_control = brail_control)
#' plot(G, method = "igraph",  weighted = TRUE)
#' 
#' @useDynLib MixedGraphs
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats runif
#' 
NULL