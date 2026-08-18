test_that("fit_lagr() supports the no-Lagrangian model", {
    x <- make_test_mcgf()
    out <- fit_lagr(x, model = "none", par_init = NULL)

    expect_equal(out$model, "none")
    expect_null(out$fit)
})

test_that("fit_lagr() fits a small one-parameter Lagrangian model", {
    x <- make_test_mcgf(80)
    x <- add_acfs(x, lag_max = 2)
    x <- add_ccfs(x, lag_max = 2)
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    fit <- fit_lagr(
        x,
        model = "lagr_tri",
        par_init = c(lambda = 0.15),
        par_fixed = c(v1 = 1.5, v2 = 0.5, k = 2)
    )

    expect_equal(fit$model, "lagr_tri")
    expect_equal(fit$par_names, "lambda")
    expect_true(is.list(fit$fit))
    expect_true(is.finite(fit$fit$objective))
    expect_true(fit$fit$par[[1]] >= 0)
    expect_true(fit$fit$par[[1]] <= 1)
})

test_that("fit_lagr() requires a previously stored base model", {
    x <- make_test_mcgf(50)
    x <- add_ccfs(x, lag_max = 2)
    attr(x, "lag") <- 2
    attr(x, "horizon") <- 1
    attr(x, "scale_time") <- 1

    expect_error(fit_lagr(
        x,
        model = "lagr_tri",
        par_init = c(lambda = 0.1),
        par_fixed = c(v1 = 1, v2 = 1, k = 2)
    ))
})
