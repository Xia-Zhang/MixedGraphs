#' glmLasso is used to fit models with lasso penalty.
#'
#' @param X is a n*p input matrix.
#' @param y is the response vector (n elements).
#' @param o is the offset vector (n elements).
#' @param lambda is the penalty weight vector, can also be a single value, we may apply the scalar to a vector in our function.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, "binomial", "gaussian" and  "poisson" are available.
#' @param support_stability is the iterations number which the support set unchanged.
#' @param thresh is the precision when the solver to stop optimizing, when the support_stability and thresh are setted at the same time, the solver would stop when any of them is achieved.
#' @param max.iter is the maximum number of iterations to be performed for the optimization.
#'
#' @return the coefficients vector
#'
#' @examples
#' X <- matrix(rnorm(500), ncol = 10)
#' y <- rbinom(50, 1, 0.6)
#' o <- rnorm(50)
#' lambda <- c(1:10)/10
#' glmLasso(X, y, o, lambda, family = "binomial", support_stability = 10, thresh = 0.5, max.iter = 1e7)
#'

glmLasso <- function(X, y, o = NULL, lambda = 1, family = c("gaussian", "binomial", "poisson"), support_stability = NULL, thresh = NULL, max.iter = 1e8) {
    list("Coef" = glmLasso_impl(X, y, o, lambda, family, support_stability, thresh, max.iter))
    if (is.atomic(lambda) && length(lambda) == 1L) {
        if (is.null(thresh)) thresh <- 1e-8
        list("Coef" = glmLasso_impl(X, y, o, lambda, family, support_stability, thresh, max.iter))
    }
    else {
        pre_coef <- NULL
        result <- sapply(lambda, function(lambda_tmp) {
            coef <- glmLasso_impl(X, y, o, lambda_tmp, family, support_stability, thresh, max.iter, pre_coef, pre_coef)
            pre_coef <<- coef
            coef
        })
        colnames(result) <- lambda
        result
    }
}

glmLasso_impl <- function(X, y, o = NULL, lambda = 1, family = c("gaussian", "binomial", "poisson"), 
    support_stability = NULL, thresh = NULL, max.iter = 1e8, init.beta = NULL, init.z = NULL, init.u = NULL) {
    X <- cbind(1, X)
    n <- nrow(X)
    p <- ncol(X)
    lambda <- append(0, lambda)
    if (is.null(o)) o <- rep(0, n)
    if (length(lambda) < p) lambda <- append(lambda, rep(lambda[2], p - length(lambda)))
    family <- tolower(family)
    family <- match.arg(family)
    if (is.null(support_stability)) {
        support_stability <- 0
    }
    else if (is.numeric(support_stability) == FALSE || support_stability <= 0) {
        stop("Invailid input support_stability!")
    }
    if (is.null(thresh)) {
        thresh <- 0.0
    }
    else if (is.numeric(thresh) == FALSE || thresh <= 0) {
        stop("Invailid input thresh!")
    }

    # warm start arguments check
    if (is.null(init.beta)) init.beta <- rep(0, p)
    else if (length(init.beta) == p - 1) init.beta <- append(0, init.beta)
    if (is.null(init.z)) init.z <- rep(0, p)
    else if (length(init.z) == p - 1) init.z <- append(0, init.z)
    if (is.null(init.u)) init.u <- rep(0, p)
    else if (length(init.u) == p - 1) init.u <- append(0, init.u)

    .Call('MixedGraphs_glmLassoCPP', PACKAGE = 'MixedGraphs', 
        X, y, o, lambda, family, support_stability, thresh, max.iter, init.beta, init.z, init.u)
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
#' glmRidge(X, y, o, lambda = 0.5, family = "binomial", thresh = 0.005, max.iter = 1e5)
#'

glmRidge <- function(X, y, o = NULL, lambda = 0.25, family = c("gaussian", "binomial", "poisson"), thresh = NULL, max.iter = 1e8) {
    if (is.atomic(lambda) && length(lambda) == 1L) {
        if (is.null(thresh)) thresh <- 1e-8
        list("Coef" = glmRidge_impl(X, y, o, lambda, family, thresh, max.iter))
    }
    else {
        pre_coef <- NULL
        if (is.null(thresh)) thresh <- 1e-1
        result <- sapply(lambda, function(lambda_tmp) {
            coef <- glmRidge_impl(X, y, o, lambda_tmp, family, thresh, max.iter, pre_coef)
            pre_coef <<- coef
            coef
        })
        colnames(result) <- lambda
        result
    }
}

glmRidge_impl <- function(X, y, o = NULL, lambda = 0.25, family = c("gaussian", "binomial", "poisson"), thresh = 1e-8, max.iter = 1e8, beta.init = NULL) {
    X <- cbind(1, X)
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(o)) o <- rep(0, n)
    if (is.null(beta.init)) beta.init <- rep(0, p)
    else if (length(beta.init) == p - 1) beta.init <- append(0, beta.init)
    family <- tolower(family)
    family <- match.arg(family)
    .Call('MixedGraphs_glmRidgeCPP', PACKAGE = 'MixedGraphs', X, y, o, beta.init, lambda, family, thresh, max.iter)
}