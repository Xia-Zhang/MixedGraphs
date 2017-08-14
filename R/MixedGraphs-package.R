#' MixedGraphs
#' 
#' A package about mixed graphical models.
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
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(c(1, 1, 0, 1), nrow = 2)
#' brail_control <- list(B = 5, tau = 0.6)
#' MixedGraph(X, crf_structure, brail_control = brail_control)
#' 
#' @useDynLib MixedGraphs
#'
#' @importFrom Rcpp evalCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats runif
#' 
NULL