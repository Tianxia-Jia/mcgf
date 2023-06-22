test_that("cor2cov() converts correlation to covariance", {
    V <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
    sd <- 1:2
    expect_equal(V, cov2cor(cor2cov(V, sd)))
})

test_that("cor2cov() errors if V is not square", {
    V <- matrix(1:6, nrow = 3)
    sd <- 1:2
    expect_error(cov2cor(cor2cov(V, sd)))
})

test_that("cor2cov() errors if lengths do not match", {
    V <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
    sd <- 1:3
    expect_error(cov2cor(cor2cov(V, sd)))
})

test_that("cor2cov() errors if V is not all non-negative", {
    V <- matrix(c(1, -0.5, 0.5, 1), ncol = 2)
    sd <- 1:3
    expect_error(cov2cor(cor2cov(V, sd)))
})
