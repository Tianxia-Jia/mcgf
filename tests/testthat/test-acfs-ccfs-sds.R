test_that("acf_rs() returns one named autocorrelation vector per regime", {
    x <- sin(seq_len(40) / 4)
    label <- factor(rep(c("A", "B"), each = 20))

    out <- acf_rs(x, label, lag_max = 3)

    expect_named(out, c("Regime A", "Regime B"))
    expect_true(all(vapply(out, length, integer(1)) == 4))
    expect_true(all(vapply(out, function(z) {
        z[[1]]
    }, numeric(1)) == 1))
    expect_equal(names(out[[1]]), paste0("lag", 0:3))
})

test_that("acfs() agrees with the mean of column-wise stats::acf", {
    x <- make_test_mcgf(50)
    got <- acfs(x, lag_max = 3)

    expected <- vapply(x, function(z) {
        stats::acf(z, lag.max = 3, plot = FALSE)$acf
    }, numeric(4))

    expect_equal(unname(got), rowMeans(expected))
    expect_equal(names(got), paste0("lag", 0:3))
})

test_that("add_acfs() stores autocorrelations and acfs<- replaces them", {
    x <- make_test_mcgf()
    x <- add_acfs(x, lag_max = 2)

    expect_equal(acfs(x), attr(x, "acfs", exact = TRUE))

    replacement <- c(
        lag0 = 1,
        lag1 = 0.25,
        lag2 = 0.1
    )
    acfs(x) <- replacement
    expect_equal(acfs(x), replacement)
})

test_that("ccfs() returns lag-zero correlations on the first slice", {
    x <- make_test_mcgf(60)
    out <- ccfs(x, lag_max = 2)

    expect_equal(dim(out), c(3, 3, 3))
    expect_equal(unname(diag(out[, , 1])), rep(1, 3), tolerance = 1e-12)
    expect_equal(out[, , 1], stats::cor(as.data.frame(x)), tolerance = 1e-12)
})

test_that("ccf_rs() returns the full negative-to-positive lag sequence", {
    x <- sin(seq_len(40) / 3)
    y <- cos(seq_len(40) / 5)
    label <- factor(rep(c("A", "B"), length.out = 40))

    out <- ccf_rs(x, y, label, lag_max = 2)

    expect_named(out, c("Regime A", "Regime B"))
    expect_true(all(vapply(out, length, integer(1)) == 5))
    expect_equal(names(out[[1]]), as.character(-2:2))
})

test_that("mcgf_rs empirical correlation summaries have regime components", {
    x <- make_test_mcgf_rs(60)

    a <- acfs(x, lag_max = 2)
    c <- ccfs(x, lag_max = 2)

    expect_named(a, c("acfs", "acfs_rs"))
    expect_named(c, c("ccfs", "ccfs_rs"))
    expect_length(a$acfs_rs, 2)
    expect_length(c$ccfs_rs, 2)
    expect_equal(dim(c$ccfs), c(3, 3, 3))
    expect_equal(dim(c$ccfs_rs[[1]]), c(3, 3, 3))
})

test_that("sds() and sd_rs() calculate standard deviations correctly", {
    x <- make_test_mcgf(40)
    expect_equal(sds(x), vapply(x, stats::sd, numeric(1)))

    xr <- make_test_mcgf_rs(40)
    label <- attr(xr, "label", exact = TRUE)
    got <- sd_rs(xr, label)

    expect_length(got, 2)
    expect_equal(got[[1]], vapply(as.data.frame(xr)[label == levels(label)[1], , drop = FALSE], stats::sd, numeric(1)))
})

test_that("add_ccfs() stores ccfs and sds together", {
    x <- make_test_mcgf(40)
    x <- add_ccfs(x, lag_max = 2)

    expect_false(is.null(attr(x, "ccfs", exact = TRUE)))
    expect_false(is.null(attr(x, "sds", exact = TRUE)))
    expect_equal(sds(x), vapply(x, stats::sd, numeric(1)))
})
