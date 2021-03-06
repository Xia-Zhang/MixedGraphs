context("Test Newton Solver")
library("glmnet")

compare_diff <- function(x1, x2, thresh = 0.5) {
    diff <- sqrt(sum((x1 - x2)^2))
    comp_result <- diff <= thresh * sqrt(length(x1))
    if (comp_result == FALSE) {
        print(x1)
        print(x2)
    }
    comp_result
}

test_that("Compare Newton Logistic solver with glmnet.", {
    n <- 20
    p <- 50
    set.seed(1)
    X <- matrix(rnorm(n * p), n, p)

    y <- rep(c(0, 1), 10)
    result_glmnet <- coef(glmnet(X, y, family = "binomial", alpha = 0), s = 0.5)
    result_ridge <- glmRidge(X, y, family = "binomial", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ridge))

    y <- rpois(n, 3)
    result_glmnet <- coef(glmnet(X, y, family = "poisson", alpha = 0), s = 0.5)
    result_ridge <- glmRidge(X, y, family = "poisson", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ridge))

    y <- rnorm(n)
    result_glmnet <- coef(glmnet(X, y, family = "gaussian", alpha = 0), s = 0.5)
    result_ridge <- glmRidge(X, y, family = "gaussian", lambda = 0.5) $ Coef
    expect_true(compare_diff(result_glmnet, result_ridge))
})

test_that("Compare Newton Logistic solver with glm.", {
    n <- 50
    p <- 10
    set.seed(1)
    X <- matrix(rnorm(n * p), n, p)

    y <- rbinom(n, 1, 0.6)
    result_glm <- coef(glm(y ~ X, family = "binomial"))
    result_ridge <- glmRidge(X, y, family = "binomial", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ridge, 1e-6))

    y <- rpois(n, 3)
    result_glm <- coef(glm(y ~ X, family = "poisson"))
    result_ridge <- glmRidge(X, y, family = "poisson", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ridge, 1e-6))

    y <- rnorm(n)
    result_glm <- coef(glm(y ~ X, family = "gaussian"))
    result_ridge <- glmRidge(X, y, family = "gaussian", lambda = 0.0) $ Coef
    expect_true(compare_diff(result_glm, result_ridge, 1e-6))
})

test_that("Test glmRidge parameter checking", {
    n <- 5
    p <- 20
    set.seed(1)
    X <- matrix(rnorm(n * p), n, p)
    y <- rnorm(n)

    expect_error(glmRidge())
    expect_error(glmRidge(X, y, family = "none"))
    expect_error(glmRidge(X, y, thresh = -1))
})

test_that("Test glmRidge with multiple lambdas input", {
    n <- 5
    p <- 20
    set.seed(1)
    X <- matrix(rnorm(n * p), n, p)
    y <- rnorm(n)

    result <- glmRidge(X, y, lambda = c(1 : 5))
    expect_equal(ncol(result), 5)
})