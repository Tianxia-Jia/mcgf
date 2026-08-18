test_that("fit_base() supports the no-base model without optimization", {
    x <- make_test_mcgf()
    out <- fit_base(x, model = "none", lag = 1)

    expect_equal(out$model, "none")
    expect_null(out$fit)
})

test_that("fit_base() fits a small one-parameter spatial model", {
    x <- make_test_mcgf(80)
    x <- add_acfs(x, lag_max = 2)
    x <- add_ccfs(x, lag_max = 2)

    fit <- fit_base(
        x,
        model = "spatial",
        lag = 1,
        par_init = c(c = 0.2),
        par_fixed = c(gamma = 0.5, nugget = 0)
    )

    expect_equal(fit$model, "spatial")
    expect_equal(fit$method, "wls")
    expect_equal(fit$par_names, "c")
    expect_true(is.list(fit$fit))
    expect_true(is.finite(fit$fit$objective))
    expect_true(fit$fit$par[[1]] >= 0)
})

test_that("fit_base() checks parameter names and lag availability", {
    x <- make_test_mcgf(50)
    x <- add_acfs(x, lag_max = 1)
    x <- add_ccfs(x, lag_max = 1)

    expect_error(
        fit_base(
            x,
            model = "spatial",
            lag = 1,
            par_init = c(bad = 1),
            par_fixed = c(gamma = 0.5, nugget = 0)
        ),
        "unknown parameters"
    )

    expect_error(
        fit_base(
            x,
            model = "spatial",
            lag = 5,
            par_init = c(c = 0.2),
            par_fixed = c(gamma = 0.5, nugget = 0)
        ),
        "recompute `acfs` and `ccfs`"
    )

    expect_error(
        fit_base(
            x,
            model = "spatial",
            method = "mle",
            lag = 1,
            par_init = c(c = 0.2),
            par_fixed = c(gamma = 0.5, nugget = 0)
        ),
        "mle is available"
    )
})
