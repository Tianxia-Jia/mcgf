test_that("mcgf() creates an mcgf data.frame with locations", {
    x <- make_test_mcgf()

    expect_s3_class(x, "mcgf")
    expect_s3_class(x, "data.frame")
    expect_false(is.mcgf_rs(x))
    expect_true(is.mcgf(x))
    expect_equal(dim(x), c(60, 3))
    expect_equal(attr(x, "locations", exact = TRUE), make_test_locations())
    expect_false(attr(x, "longlat", exact = TRUE))
})

test_that("mcgf() accepts precomputed distances", {
    dat <- make_test_data(20)
    d <- make_test_dists()

    x <- suppressMessages(mcgf(dat, dists = d))

    expect_s3_class(x, "mcgf")
    expect_equal(dists(x)$h, d$h)
    expect_equal(dists(x)$h1, d$h1)
    expect_equal(dists(x)$h2, d$h2)
})

test_that("mcgf() validates its core inputs", {
    dat <- make_test_data(10)
    loc <- make_test_locations()
    d <- make_test_dists()

    expect_error(mcgf(as.character(dat$S1), locations = loc))
    expect_error(mcgf(transform(dat, S1 = NA_real_), locations = loc))
    expect_error(mcgf(transform(dat, S1 = as.character(S1)), locations = loc))
    expect_error(mcgf(dat))
    expect_error(mcgf(dat, locations = loc, dists = d))
    expect_error(mcgf(dat, locations = loc[-1, , drop = FALSE]))
    expect_error(mcgf(dat, locations = loc, time = c(1:9, 11)))
    expect_error(mcgf(dat, locations = loc, time = 10:1))
})

test_that("mcgf() preserves an equally spaced supplied time index", {
    dat <- make_test_data(8)
    time <- seq.Date(as.Date("2026-01-01"),
        by = "day",
        length.out = 8
    )

    x <- mcgf(dat,
        locations = make_test_locations(),
        time = time,
        longlat = FALSE
    )

    expect_equal(rownames(x), as.character(time))
})
