#' Estimation function with lasso penalty
#'
#' @description glmLasso is used to fit models with lasso penalty.
#'
#' @param X is a n*p data matrix.
#' @param y is the response vector of length n.
#' @param o is the offset vector of length n.
#' @param lambda is the penalty parameter, can be either a scalar or a vector. If a vector, would apply each element to do a regression.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, the character string can be "binomial", "gaussian" or "poisson".
#' @param support_stability decides the stop criterian. When the support set stay unchanged for support_stability iterations, it assume that the function is stable enough to stop.
#' @param thresh is the threshold of the diffrence between adjacent iteration coefficents. Parameter thresh and support_stability can not be set at the same time.
#' @param max.iter is the maximum number of iterations to be performed for the optimization routine.
#' @param intercept is a boolean value, which indicates whether intercept will be fitted in the function. Default value is TRUE.
#'
#' @return the coefficients vector
#'
#' @examples
#' n <- 5
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rbinom(n, 1, 0.6)
#' o <- rnorm(n)
#' lambda <- c(1:5)/10
#' glmLasso(X, y, o, lambda, family = "binomial", support_stability = 10, 
#'          max.iter = 1e7, intercept = TRUE)
#'
#' @export

glmLasso <- function(X, y, o = NULL, lambda = 1, family = c("gaussian", "binomial", "poisson"),
                     support_stability = NULL, thresh = NULL, max.iter = 1e8, intercept = TRUE) {
    names <- colnames(X, do.NULL = FALSE, prefix = "v")
    if (intercept) names <- c("(Intercept)", names)
    if (is.atomic(lambda) && length(lambda) == 1L) {
        if (is.null(thresh) && is.null(support_stability)) thresh <- 1e-8
        coefficients = glmLasso_impl(X, y, o, weight = lambda, NULL, family, support_stability, thresh, max.iter, intercept)
        rownames(coefficients) <- names
        list("Coef" = coefficients)
    }
    else {
        pre_coef <- NULL
        result <- glmLasso_impl(X, y, o, NULL, lambda, family, support_stability, thresh, max.iter, intercept, pre_coef, pre_coef)
        colnames(result) <- lambda
        rownames(result) <- names
        result
    }
}

glmLasso_impl <- function(X, y, o = NULL, weight = NULL, lambda = NULL, family = c("gaussian", "binomial", "poisson"), support_stability = NULL,
                          thresh = NULL, max.iter = 1e8, intercept = TRUE, init.beta = NULL, init.z = NULL, init.u = NULL) {
    if (intercept) {
        X <- cbind(1, X)
    }
    n <- nrow(X)
    p <- ncol(X)

    if (is.null(o)) o <- rep(0, n)
    if (is.null(lambda) && is.null(weight)) {
        stop("At least one of weight and lambda should be exist!")
    }
    if (is.null(lambda)) {
        if (intercept) {
            weight <- append(0, weight)
            weight <- append(weight, rep(weight[2], p - length(weight)))
        }
        else {
            weight <- rep(weight, p)
        }
        lambda <- numeric()
    }
    else {
        weight <- numeric()
    }
    family <- tolower(family)
    family <- match.arg(family)

    if (!is.null(support_stability) && !is.null(thresh)) {
        stop("support_stability and thresh cannot be set at the same time!")
    }
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

    glmLassoCPP(X, y, o, weight, lambda, family, support_stability, thresh, max.iter, init.beta, init.z, init.u, intercept)
}

#' Estimation function with ridge penalty
#'
#' @description glmRidge is used to fit models with ridge penalty.
#'
#' @param X is a n*p data matrix.
#' @param y is the response vector of length n.
#' @param o is the offset vector of length n.
#' @param lambda is the penalty parameter, can be either a scalar or a vector. If a vector, would apply each element to do a regression.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, the character string can be "binomial", "gaussian" or "poisson".
#' @param thresh is the threshold to stop the solver. The comparison would perform between it and the difference of the the adjacent iteration coefficents Euclidean norm .
#' @param max.iter is the maximum number of iterations to be performed for the optimization.
#' @param intercept is a boolean value, which indicate whether intercept will be fitted in the function. Default value is TRUE.
#'
#' @return the coefficients vector
#'
#' @examples
#' n <- 50
#' p <- 10
#' X <- matrix(rnorm(n * p), n, p)
#' y <- rbinom(n, 1, 0.6)
#' o <- rnorm(n)
#' glmRidge(X, y, o, lambda = 0.5, family = "binomial", thresh = 0.005, 
#'          max.iter = 1e5, intercept = TRUE)
#'
#' @export

glmRidge <- function(X, y, o = NULL, lambda = 0.25, family = c("gaussian", "binomial", "poisson"),
                     thresh = NULL, max.iter = 1e8, intercept = TRUE) {
    names <- colnames(X, do.NULL = FALSE, prefix = "v")
    if (intercept) names <- c("(Intercept)", names)
    if (is.atomic(lambda) && length(lambda) == 1L) {
        if (is.null(thresh)) thresh <- 1e-8
        coefficients <- glmRidge_impl(X, y, o, lambda, family, thresh, max.iter, intercept)
        rownames(coefficients) <- names
        list("Coef" = coefficients)
    }
    else {
        pre_coef <- NULL
        if (is.null(thresh)) thresh <- 1e-1
        result <- glmRidge_impl(X, y, o, lambda, family, thresh, max.iter, intercept, pre_coef)
        colnames(result) <- lambda
        rownames(result) <- names
        result
    }
}

glmRidge_impl <- function(X, y, o = NULL, lambda = 0.25, family = c("gaussian", "binomial", "poisson"),
                          thresh = 1e-8, max.iter = 1e8, intercept = TRUE, beta.init = NULL) {
    intercept_mean <- 0
    family <- tolower(family)
    family <- match.arg(family)
    if (intercept) {
        X <- cbind(1, X)
        if (family == "gaussian") {
            intercept_mean <- mean(y)
            y <- y -  mean(y)
        }
    }
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(o)) o <- rep(0, n)
    if (is.null(beta.init)) beta.init <- rep(0, p)
    else if (length(beta.init) == p - 1) beta.init <- append(0, beta.init)
    if (is.numeric(thresh) == FALSE || thresh <= 0) {
        stop("Invailid input thresh!")
    }
    result <- glmRidgeCPP(X, y, o, beta.init, as.vector(lambda), family, thresh, max.iter, intercept)
    if (intercept && family == "gaussian") {
        result[1, ] <- result[1, ] + intercept_mean
    }
    result
}