#' glmLasso is used to fit models with lasso penalty.
#'
#' @param X is a n*p input matrix.
#' @param y is the response vector (n elements).
#' @param o is the offset vector (n elements).
#' @param lambda is the penalty weight vector, can also be a single value, we may apply the scalar to a vector in our function.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, "binomial", "gaussian" and  "poisson" are available.
#' @param KLB is the iterations number which the support set unchanged.
#' @param thresh is the precision when the solver to stop optimizing, when the KLB and thresh are setted at the same time, the solver would stop when any of them is achieved.
#' @param max.iter is the maximum number of iterations to be performed for the optimization.
#'
#' @return the coefficients vector
#'
#' @examples
#' X <- matrix(rnorm(500), ncol = 10)
#' y <- rbinom(50, 1, 0.6)
#' o <- rnorm(50)
#' lambda <- c(1:10)/10
#' glmLasso(X, y, o, lambda, family = "binomial", KLB = 10, thresh = 0.5, max.iter = 1e7)
#'

glmLasso <- function(X, y, o = NULL, lambda = 1, family = "gaussian", KLB = NULL, thresh = NULL, max.iter = 1e8) {
	X <- cbind(1, X)
	n <- nrow(X)
	p <- ncol(X)
	lambda <- append(0, lambda)
	if (is.null(o)) o <- rep(0, n)
	if (length(lambda) < p) lambda <- append(lambda, rep(lambda[2], p - length(lambda)))
	family <- tolower(family)
	if (family == "binomial") 
		family <- "logistic"
	else if (is.element(family, c('gaussian', 'logistic', 'poisson')) == FALSE) 
		stop("The family should be in c('gaussian', 'logistic', 'poisson').")
	if (is.null(KLB)) {
		KLB <- 0
	}
	else if (is.numeric(KLB) == FALSE || KLB <= 0) {
		stop("Invailid input KLB!")
	}
	if (is.null(thresh)) {
		thresh <- 0.0
	}
	else if (is.numeric(thresh) == FALSE || thresh <= 0) {
		stop("Invailid input thresh!")
	}
	.Call('MixedGraphs_glmLassoCPP', PACKAGE = 'MixedGraphs', X, y, o, lambda, family, KLB, thresh, max.iter)
}


#' glmRidge is used to fit models with ridge penalty.
#'
#' @param X is a n*p input matrix.
#' @param y is the response vector (n elements).
#' @param o is the offset vector (n elements).
#' @param lambda is a single value of ridge penalty.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, "binomial", "gaussian" and  "poisson" are available.
#' @param thresh indicates when to stop the solver. glmRidge will check the Euclidean norm of the change.
#' @param max.iter is the maximum number of iterations to be performed for the optimization.
#'
#' @return the coefficients vector
#'
#' @examples
#' X <- matrix(rnorm(500), ncol = 10)
#' y <- rbinom(50, 1, 0.6)
#' o <- rnorm(50)
#'
#' glmRidge(X, y, o, lambda = 0.5, family = "binomial", thresh = 0.005, max.iter = 1e5)
#'

glmRidge <- function(X, y, o = NULL, lambda = 0.25, family = "gaussian", thresh = 1e-8, max.iter = 1e8) {
	X <- cbind(1, X)
	n <- nrow(X)
	p <- ncol(X)
	if (is.null(o)) o <- rep(0, n)
	if (family == "binomial")
		family <- "logistic"
	else if (is.element(family, c('gaussian', 'logistic', 'poisson')) == FALSE)
		stop("The family should be in c('gaussian', 'logistic', 'poisson').")

	.Call('MixedGraphs_glmRidgeCPP', PACKAGE = 'MixedGraphs', X, y, o, lambda, family, thresh, max.iter)
}