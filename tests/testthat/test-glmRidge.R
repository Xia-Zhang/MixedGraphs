context("Test Newton Solver")
library("glmnet")

compare_diff <- function(x1, x2, thresh = 0.5) {
    diff <- sqrt(sum((x1 - x2)^2))
    diff <= thresh * sqrt(length(x1))
}

test_that("Compare Newton Logistic solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "binomial", alpha = 0), s = 0.5)
    result_ADMM <- glmRidge(X, y, family = "binomial", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ADMM))
})

test_that("Compare Newton Poisson solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rpois(50, 3)

    result_glmnet <- coef(glmnet(X, y, family = "poisson", alpha = 0), s = 0.5)
    result_ADMM <- glmRidge(X, y, family = "poisson", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ADMM))
})

test_that("Compare Newton Gaussian solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rnorm(50)

    result_glmnet <- coef(glmnet(X, y, family = "gaussian", alpha = 0), s = 0.5)
    result_ADMM <- glmRidge(X, y, family = "gaussian", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ADMM))
})

test_that("Compare Newton Logistic solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glm <- coef(glm(y ~ X, family = "binomial"))
    result_ADMM <- glmRidge(X, y, family = "binomial", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM, 1e-6))
})

test_that("Compare Newton Poisson solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rpois(50, 3)

    result_glm <- coef(glm(y ~ X, family = "poisson"))
    result_ADMM <- glmRidge(X, y, family = "poisson", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM, 1e-6))
})

test_that("Compare Newton Gaussian solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rnorm(50)

    result_glm <- coef(glm(y ~ X, family = "gaussian"))
    result_ADMM <- glmRidge(X, y, family = "gaussian", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM, 1e-6))
})