test_that("ccov() constructs an empirical joint covariance matrix", {
    x <- make_test_mcgf(60)
    x <- add_ccfs(x, lag_max = 2)
    attr(x, "lag") <- 2
    attr(x, "horizon") <- 1

    out <- ccov(x, model = "empirical")

    expect_equal(dim(out), c(9, 9))
    expect_true(isSymmetric(out, tol = 1e-10))
})

test_that("ccov() constructs a base-model joint correlation matrix", {
    x <- make_test_mcgf(60)
    sds(x) <- sds(x)
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    out <- ccov(x, model = "base", cor = TRUE)

    expect_equal(dim(out), c(9, 9))
    expect_equal(unname(diag(out)[1:3]), rep(1, 3))
})

test_that("krige() forecasts from empirical covariance", {
    x <- make_test_mcgf(60)
    x <- add_ccfs(x, lag_max = 2)

    pred <- krige(x, model = "empirical", lag = 2, horizon = 1)

    expect_equal(dim(pred), c(60, 3))
    expect_true(any(is.finite(pred)))
})

test_that("krige() returns prediction intervals when requested", {
    x <- make_test_mcgf(60)
    x <- add_ccfs(x, lag_max = 2)

    pred <- krige(x, model = "empirical", lag = 2, horizon = 1, interval = TRUE)

    expect_named(pred, c("fit", "lower", "upper"))
    expect_equal(dim(pred$fit), c(60, 3))
    expect_true(all(pred$lower <= pred$upper, na.rm = TRUE))
})

test_that("krige() validates required fitted models", {
    x <- make_test_mcgf()

    expect_error(krige(x, model = "base"), "base model missing")
    expect_error(krige(x, model = "all"), "Lagrangian model missing")
})
