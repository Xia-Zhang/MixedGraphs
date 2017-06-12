context("Test ADMM Solver")


compare_supp <- function(x1, x2) {
    count <- 0
    for (i in 1:length(x1)) {
        if (x1[i] == 0 && x2[i] == 0 || x1[i] != 0 && x2[i] != 0) {
            count <- count + 1
        }
    }
    length(x1) == count
}

test_that("Compare ADMM Logistic solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)
    library("glmnet")

    result_glmnet <- coef(glmnet(X, y, family = "binomial"), s = 1)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "logistic", 1) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})


test_that("Compare ADMM Poisson solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rpois(50, 3)

    result_glmnet <- coef(glmnet(X, y, family = "poisson"), s = 1)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "poisson", 1) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})

test_that("Compare ADMM Gaussian solver with glmnet.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rnorm(50)

    result_glmnet <- coef(glmnet(X, y, family = "gaussian"), s = 1)
    result_ADMM <- .Call('MixedGraphs_testADMM', PACKAGE = 'MixedGraphs', cbind(1, X), y, "gaussian", 1) $ Result
    expect_true(compare_supp(result_glmnet, result_ADMM))
})