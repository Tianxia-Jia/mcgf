test_that("cor2cov() applies marginal standard deviations", {
    cor <- matrix(c(1, 0.25, 0.25, 1), 2)
    sd <- c(2, 3)

    expected <- diag(sd) %*% cor %*% diag(sd)
    expect_equal(cor2cov(cor, sd), expected)
})

test_that("cor2cov_ar() applies cor2cov() to every lag slice", {
    cor <- array(c(1, 0.25, 0.25, 1, 1, 0.10, 0.10, 1), dim = c(2, 2, 2))
    sd <- c(2, 3)

    got <- cor2cov_ar(cor, sd)
    expect_equal(got[, , 1], cor2cov(cor[, , 1], sd))
    expect_equal(got[, , 2], cor2cov(cor[, , 2], sd))
})

test_that("cov_joint() creates the expected block Toeplitz matrix", {
    c0 <- matrix(c(1, 0.2, 0.2, 1), 2)
    c1 <- matrix(c(0.5, 0.1, 0.3, 0.4), 2)
    cov_ar <- array(c(c0, c1), dim = c(2, 2, 2))

    got <- mcgf:::cov_joint(cov_ar)
    expected <- rbind(cbind(c0, c1), cbind(t(c1), c0))

    expect_equal(got, expected)
})

test_that("cov_par() returns conditional weights and covariance", {
    joint <- matrix(c(1, 0.4, 0.4, 1), 2)

    out <- mcgf:::cov_par(joint, n_var = 1, joint = TRUE)

    expect_equal(out$weights, matrix(0.4, 1, 1))
    expect_equal(out$cov_curr, matrix(1 - 0.4^2, 1, 1))
})

test_that("to_ar() repeats a distance matrix and constructs lag arrays", {
    h <- matrix(c(0, 1, 1, 0), 2)
    out <- mcgf:::to_ar(h, lag_max = 2)

    expect_equal(dim(out$h_ar), c(2, 2, 3))
    expect_equal(out$h_ar[, , 1], h)
    expect_equal(out$h_ar[, , 3], h)
    expect_equal(unique(c(out$u_ar[, , 1])), 0)
    expect_equal(unique(c(out$u_ar[, , 2])), 1)
    expect_equal(unique(c(out$u_ar[, , 3])), 2)
})

test_that("mat_inv() returns a generalized inverse when necessary", {
    x <- matrix(c(1, 1, 1, 1), 2)
    inv <- mcgf:::mat_inv(x)

    expect_equal(x %*% inv %*% x, x, tolerance = 1e-8)
})
