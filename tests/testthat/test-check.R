# Check is_numeric_scalar()

test_that("is_numeric_scalar() returns FALSE if x is not a scalar", {
    expect_equal(FALSE, is_numeric_scalar(x = 1:2))
    expect_equal(FALSE, is_numeric_scalar(x = "a"))
    expect_equal(FALSE, is_numeric_scalar(x = NaN))
    expect_equal(FALSE, is_numeric_scalar(x = NA))
    expect_equal(FALSE, is_numeric_scalar(x = NULL))
})

# Check check_dist()

test_that("check_dist() errors if any element in x is NA or NaN", {
    expect_error(check_dist(x = NA))

    h <- matrix(c(NA, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(check_dist(x = h))

    h <- array(matrix(c(NaN, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3),
        dim = c(3, 3, 3)
    )
    expect_error(check_dist(x = h))
})

test_that("check_dist() errors if any element in x is negative", {
    expect_error(check_dist(x = -1))

    h <- matrix(c(-1, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(check_dist(x = h))

    h <- array(matrix(c(-0.1, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3),
        dim = c(3, 3, 3)
    )
    expect_error(check_dist(x = h))
})

test_that("check_dist() errors if x is not a matrix or 3d array", {
    h <- matrix(c(0, 1, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(check_dist(x = h))

    h <- array(matrix(c(0, 1, 0, 5, 0, 1, 0, 1, 0), nrow = 3),
        dim = c(3, 3, 3)
    )
    expect_error(check_dist(x = h))
})
