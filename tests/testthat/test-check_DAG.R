context("Test DAG check function")

test_that("make sure the DAG judge result", {
    crf_structure = matrix(0, 10, 10)
    expect_equal(check_DAG(crf_structure), TRUE)

    crf_structure = matrix(1, 10, 10)
    expect_equal(check_DAG(crf_structure), FALSE)

    crf_structure = matrix(c(0, 1, 0, 0), 2, 2)
    expect_equal(check_DAG(crf_structure), TRUE)

    crf_structure = matrix(c(0, 1, 1, 0), 2, 2)
    expect_equal(check_DAG(crf_structure), FALSE)

    crf_structure = matrix(c(1, 0, 1, 
                             1, 1, 1, 
                             0, 0, 1), 3, 3)
    expect_equal(check_DAG(crf_structure), TRUE)

    crf_structure = matrix(c(1, 1, 1, 
                             1, 1, 1, 
                             0, 0, 1), 3, 3)
    expect_equal(check_DAG(crf_structure), FALSE)

    crf_structure = matrix(c(1, 1, 0, 
                             0, 1, 1, 
                             1, 0, 1), 3, 3) 
    expect_equal(check_DAG(crf_structure), FALSE)
})
