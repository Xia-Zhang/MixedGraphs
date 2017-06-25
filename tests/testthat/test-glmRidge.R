context("Test ADMM Solver")
library("glmnet")

compare_supp <- function(x1, x2) {
    comp_result <- all(sign(x1) == sign(x2))
    if (comp_result == FALSE) {
        print(x1)
        print(x2)
    }
    comp_result
}

compare_diff <- function(x1, x2, thresh = 0.5) {
    diff <- sqrt(sum((x1 - x2)^2))
    diff <= thresh * sqrt(length(x1))
}

test_that("Compare support set of glmLasso and glmnet, when set KLB to 5 and family to binomial.", {
    X <- matrix(rnorm(1000), ncol = 10)
    y <- rbinom(100, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "binomial", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- glmLasso(X, y, family = "binomial", lambda = 0.5) $ Coef
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare support set of glmLasso and glmnet, when set KLB to 5 and family to poisson.", {
    X <- matrix(rnorm(1000), ncol = 5)
    y <- rpois(200, 3)
    result_glmnet <- coef(glmnet(X, y, family = "poisson", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- glmLasso(X, y, family = "poisson", lambda = 0.5) $ Coef
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare support set of glmLasso and glmnet, when set KLB to 5 and family to gaussian.", {
    X <- matrix(rnorm(1000), ncol = 10)
    y <- rnorm(100)

    result_glmnet <- coef(glmnet(X, y, family = "gaussian", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- glmLasso(X, y, family = "gaussian", lambda = 0.5) $ Coef
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare the glmLasso and glm, when set thresh to 0.1 and family to binomial.", {
    X <- matrix(rnorm(250), ncol = 5)
    y <- rbinom(50, 1, 0.6)

    result_glm <- coef(glm(y ~ X, family = "binomial"))
    result_ADMM <- glmLasso(X, y, family = "binomial", lambda = 0, thresh = 1e-2) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM))
})

test_that("Compare the glmLasso and glm, when set thresh to 0.1 and family to poisson.", {
    X <- matrix(rnorm(250), ncol = 5)
    y <- rpois(50, 3)

    result_glm <- coef(glm(y ~ X, family = "poisson"))
    result_ADMM <- glmLasso(X, y, family = "poisson", lambda = 0, thresh = 1e-2) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM))
})

test_that("Compare the glmLasso and glm, when set thresh to 0.1 and family to gaussian.", {
    X <- matrix(rnorm(250), ncol = 5)
    y <- rnorm(50)

    result_glm <- coef(glm(y ~ X, family = "gaussian"))
    result_ADMM <- glmLasso(X, y, family = "gaussian", lambda = 0, thresh = 1e-2) $ Coef
    expect_true(compare_diff(result_glm, result_ADMM))
})