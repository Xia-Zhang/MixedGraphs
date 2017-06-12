context("Test Newton Solver")


norm_vec <- function(x) sqrt(sum(x^2))

compare_diff <- function(x1, x2) {
    diff <- norm_vec(x1 - x2)
    epsilon <- 0.5
    diff <= epsilon * sqrt(length(x1))
}

test_that("Compare Newton Logistic solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "binomial", alpha = 0), s = 0.5)
    result_ADMM <- testNewton(cbind(1, X), y, "logistic", 0.5) $ Result
    expect_true(compare_diff(result_glmnet, result_ADMM))
})

test_that("Compare Newton Poisson solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "poisson", alpha = 0), s = 0.5)
    result_ADMM <- testNewton(cbind(1, X), y, "poisson", 0.5) $ Result
    expect_true(compare_diff(result_glmnet, result_ADMM))
})

test_that("Compare Newton Gaussian solver with glm.", {
    X <- matrix(rnorm(500), ncol = 10)
    y <- rbinom(50, 1, 0.6)

    result_glmnet <- coef(glmnet(X, y, family = "gaussian", alpha = 0), s = 0.5)
    result_ADMM <- testNewton(cbind(1, X), y, "Gaussian", 0.5) $ Result
    expect_true(compare_diff(result_glmnet, result_ADMM))
})