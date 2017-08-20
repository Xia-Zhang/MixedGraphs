context("Test family guess function")

test_that("check if family guess function work", {
    n <- 3
    p <- 4

    X1 <- matrix(rnorm(n * p), n, p)
    expect_equal(guess_family(X1), "gaussian")

    X2 <- matrix(rbinom(n * p, 1, 0.6), n, p)
    expect_equal(guess_family(X2), "binomial")

    X3 <- matrix(rpois(n * p, 3), n, p)
    expect_equal(guess_family(X3), "poisson")

    expect_equal(guess_family(list(X3, X1, X2)), c("poisson", "gaussian", "binomial"))
})
