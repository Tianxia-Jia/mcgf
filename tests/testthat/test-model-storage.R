test_that("add_base() stores a fitted base model and its covariance", {
    x <- make_test_mcgf()
    fit <- make_mock_base_fit(x, lag = 2)

    y <- add_base(x, fit_base = fit)

    expect_equal(attr(y, "base", exact = TRUE), "sep")
    expect_equal(attr(y, "lag", exact = TRUE), 2)
    expect_equal(attr(y, "horizon", exact = TRUE), 1)
    expect_equal(dim(attr(y, "base_res", exact = TRUE)$cor_base), c(3, 3, 3))
})

test_that("base<- stores a complete user-supplied base model", {
    x <- make_test_mcgf()
    fit <- make_mock_base_fit(x, lag = 2)

    suppressMessages(base(x) <- fit)

    expect_equal(attr(x, "base", exact = TRUE), "sep")
    expect_equal(attr(x, "base_res", exact = TRUE)$cor_base, fit$cor_base)
    expect_equal(attr(x, "fit_base_raw", exact = TRUE), fit)
})

test_that("add_lagr() stores the Lagrangian covariance after a base model", {
    x <- make_test_mcgf()
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))
    fit <- make_mock_lagr_fit(x)

    y <- add_lagr(x, fit_lagr = fit)

    expect_equal(attr(y, "lagr", exact = TRUE), "lagr_tri")
    expect_equal(dim(attr(y, "lagr_res", exact = TRUE)$cor_lagr), c(3, 3, 3))
})

test_that("lagr<- stores a user-supplied Lagrangian model", {
    x <- make_test_mcgf()
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))
    fit <- make_mock_lagr_fit(x)

    suppressMessages(lagr(x) <- fit)

    expect_equal(attr(x, "lagr", exact = TRUE), "lagr_tri")
    expect_equal(attr(x, "lagr_res", exact = TRUE)$cor_lagr, fit$cor_lagr)
})

test_that("model() and print.mcgf() provide readable summaries", {
    x <- make_test_mcgf()
    x <- add_base(x, fit_base = make_mock_base_fit(x, lag = 2))

    out_model <- capture.output(model(x, model = "base"))
    out_print <- capture.output(print(x, "dists"))

    expect_true(any(grepl("Base model", out_model, fixed = TRUE)))
    expect_true(length(out_print) > 0)
})
