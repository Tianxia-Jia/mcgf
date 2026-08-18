test_that("mcgf_rs() creates the expected inherited classes and labels", {
    x <- make_test_mcgf_rs()

    expect_s3_class(x, "mcgf_rs")
    expect_true(is.mcgf_rs(x))
    expect_true(is.mcgf(x))
    expect_equal(levels(attr(x, "label", exact = TRUE)), c("A", "B"))
    expect_length(attr(x, "label", exact = TRUE), nrow(x))
})

test_that("mcgf_rs() requires one label per observation", {
    dat <- make_test_data(20)

    expect_error(suppressMessages(
        mcgf_rs(
            dat,
            locations = make_test_locations(),
            label = rep(c("A", "B"), 9),
            longlat = FALSE
        )
    ), "length of `label`")
})

test_that("as.mcgf_rs() converts an mcgf object without changing the data", {
    x <- make_test_mcgf(20)
    label <- rep(c(1, 2), length.out = nrow(x))

    y <- as.mcgf_rs(x, label = label)

    expect_true(is.mcgf_rs(y))
    expect_equal(c(as.data.frame(y)), c(as.data.frame(x)))
    expect_equal(
        as.character(attr(y, "label", exact = TRUE)),
        as.character(factor(label))
    )
})
