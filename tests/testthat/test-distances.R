test_that("find_dists() returns consistent Euclidean signed distances", {
    loc <- make_test_locations()
    d <- find_dists(loc, longlat = FALSE)

    expect_named(d, c("h", "h1", "h2"))
    expect_equal(dim(d$h), c(3, 3))
    expect_true(isSymmetric(d$h))
    expect_equal(d$h, sqrt(d$h1^2 + d$h2^2))
    expect_equal(d$h1, -t(d$h1))
    expect_equal(d$h2, -t(d$h2))
    expect_equal(unname(diag(d$h)), rep(0, 3))
    expect_equal(rownames(d$h), rownames(loc))
})

test_that("find_dists_new() appends new locations to the same distance system", {
    loc <- make_test_locations()
    loc_new <- matrix(c(1, 0.75), nrow = 1, dimnames = list("S4", c("x", "y")))

    d_old <- find_dists(loc, longlat = FALSE)
    d_new <- find_dists_new(loc, loc_new, longlat = FALSE)

    expect_equal(dim(d_new$h), c(4, 4))
    expect_equal(d_new$h[1:3, 1:3], d_old$h)
    expect_equal(d_new$h1[1:3, 1:3], d_old$h1)
    expect_equal(d_new$h2[1:3, 1:3], d_old$h2)
    expect_equal(rownames(d_new$h), c("S1", "S2", "S3", "S4"))
})

test_that("find_dists() and find_dists_new() reject malformed coordinates", {
    expect_error(find_dists(matrix(1:9, ncol = 3)))
    expect_error(find_dists_new(make_test_locations(), c(1, 2, 3)))

    loc_new <- rbind(S1 = c(1, 1), S4 = c(2, 2))
    expect_error(
        find_dists_new(make_test_locations(), loc_new, longlat = FALSE),
        "duplicate row names"
    )
})

test_that("rdists() is reproducible and respects names and scale", {
    set.seed(11)
    d1 <- rdists(4, names = LETTERS[1:4], scale = 2)
    set.seed(11)
    d2 <- rdists(4, names = LETTERS[1:4], scale = 2)

    expect_equal(d1, d2)
    expect_equal(dim(d1$h), c(4, 4))
    expect_equal(rownames(d1$h), LETTERS[1:4])
    expect_error(rdists(4, names = LETTERS[1:3]))
})

test_that("dists() calculates and replacement stores validated distances", {
    x <- make_test_mcgf()
    d <- dists(x)

    expect_equal(d$h, make_test_dists()$h)

    x2 <- x
    dists(x2) <- d
    expect_equal(attr(x2, "dists", exact = TRUE)$h, d$h)
})
