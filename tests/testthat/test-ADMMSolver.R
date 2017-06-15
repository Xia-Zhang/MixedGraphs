context("Test ADMM Solver")
library("glmnet")

compare_supp <- function(x1, x2) {
    all(sign(x1) == sign(x2))
}

test_that("Compare ADMM Logistic solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "binomial", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "logistic", 0.5) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare ADMM Poisson solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rpois(50, 3)

    result_glmnet <- coef(glmnet(X, y, family = "poisson", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "poisson", 0.5) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare ADMM Gaussian solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rnorm(50)

    result_glmnet <- coef(glmnet(X, y, family = "gaussian", standardize=FALSE, standardize.response=FALSE), s = 0.5)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "gaussian", 0.5) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})