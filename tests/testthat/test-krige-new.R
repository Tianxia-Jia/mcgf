test_that("krige_new() forecasts at a new location using a stored base fit", {
    x <- make_test_mcgf(60)
    sds(x) <- sds(x)
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    pred <- krige_new(x, locations_new = c(1, 0.75), model = "base")

    expect_equal(dim(pred), c(60, 4))
    expect_equal(dimnames(pred)[[2]][1:3], colnames(x))
    expect_true(any(is.finite(pred[, 4])))
})

test_that("krige_new() can return intervals at new locations", {
    x <- make_test_mcgf(60)
    sds(x) <- sds(x)
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    pred <- krige_new(
        x,
        locations_new = c(1, 0.75),
        model = "base",
        interval = TRUE
    )

    expect_named(pred, c("fit", "lower", "upper"))
    expect_equal(dim(pred$fit), c(60, 4))
    expect_true(all(pred$lower <= pred$upper, na.rm = TRUE))
})

test_that("krige_new() validates location and distance inputs", {
    x <- make_test_mcgf(30)
    sds(x) <- sds(x)
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    expect_error(krige_new(x, model = "base"), "must provide either")

    expect_error(
        krige_new(
            x,
            locations_new = c(1, 1),
            dists_new = find_dists_new(make_test_locations(),
                c(1, 1),
                longlat = FALSE
            ),
            model = "base"
        ),
        "do not provide both"
    )
})
