context("Test ADMM Solver")
library("glmnet")

compare_supp <- function(x1, x2) {
    comp_result <- all(abs(sign(x1)) == abs(sign(x2)))
    if (comp_result == FALSE) {
        print(x1)
        print(x2)
    }
    comp_result
}

compare_diff <- function(x1, x2, thresh = 0.5) {
    diff <- sqrt(sum((x1 - x2)^2))
    comp_result <- diff <= thresh * sqrt(length(x1))
    if (comp_result == FALSE) {
        print(x1)
        print(x2)
    }
    comp_result
}

test_that("Compare support set of glmLasso and glmnet", {
    n <- 20
    p <- 5
    X <- matrix(rnorm(n * p), n, p)

    y <- rep(c(0, 1), 10)
    result_glmnet <- coef(glmnet(X, y, family = "binomial", standardize=FALSE, standardize.response=FALSE), s = 1)
    result_lasso <- glmLasso(X, y, family = "binomial", lambda = 1, support_stability = 10) $ Coef
    expect_true(compare_supp(result_glmnet, result_lasso))

    y <- rpois(n, 3)
    result_glmnet <- coef(glmnet(X, y, family = "poisson", standardize=FALSE, standardize.response=FALSE), s = 1)
    result_lasso <- glmLasso(X, y, family = "poisson", lambda = 1, support_stability = 20) $ Coef
    expect_true(compare_supp(result_glmnet, result_lasso))

    y <- rnorm(n)
    result_glmnet <- coef(glmnet(X, y, family = "gaussian", standardize=FALSE, standardize.response=FALSE), s = 1)
    result_lasso <- glmLasso(X, y, family = "gaussian", lambda = 1, support_stability = 20) $ Coef
    expect_true(compare_supp(result_glmnet, result_lasso))
})

test_that("Compare the glmLasso and glm", {
    n <- 200
    p <- 5
    X <- matrix(rnorm(n * p), n, p)

    y <- rbinom(n, 1, 0.6)
    result_glm <- coef(glm(y ~ X, family = "binomial"))
    result_lasso <- glmLasso(X, y, family = "binomial", lambda = 0, thresh = 1e-6) $ Coef
    expect_true(compare_diff(result_glm, result_lasso))

    y <- rpois(n, 3)
    result_glm <- coef(glm(y ~ X, family = "poisson"))
    result_lasso <- glmLasso(X, y, family = "poisson", lambda = 0, thresh = 1e-6) $ Coef
    expect_true(compare_diff(result_glm, result_lasso))

    y <- rnorm(n)
    result_glm <- coef(glm(y ~ X, family = "gaussian"))
    result_lasso <- glmLasso(X, y, family = "gaussian", lambda = 0, thresh = 1e-6) $ Coef
    expect_true(compare_diff(result_glm, result_lasso))
})

test_that("Test glmLasso parameter checking", {
    n <- 5
    p <- 20
    X <- matrix(rnorm(n * p), n, p)
    y <- rnorm(n)

    expect_error(glmLasso())
    expect_error(glmLasso(X, y, family = "none"))
    expect_error(glmLasso(X, y, support_stability = 20, thresh = 1e-6))
    expect_error(glmLasso(X, y, support_stability = -2))
    expect_error(glmLasso(X, y, thresh = "thresh"))
})

test_that("Test glmLasso with multiple lambdas input", {
    n <- 5
    p <- 20
    X <- matrix(rnorm(n * p), n, p)
    y <- rnorm(n)

    result <- glmLasso(X, y, lambda = c(1 : 5))
    expect_equal(ncol(result), 5)
})