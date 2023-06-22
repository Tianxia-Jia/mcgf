test_that("cor_cauchy() calculates Cauchy correlation", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    corr <- cor_cauchy(x, a = 5, alpha = 0.5, nu = 3)
    expect_equal(corr, .cor_cauchy(x = x, a = 5, alpha = 0.5, nu = 3))
})

test_that("cor_cauchy() errors if nugget > 0 but is.dist = FALSE", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(cor_cauchy(x = x, a = 5, alpha = 0.5, nugget = 0.5,
                            is.dist = FALSE))
})

test_that("cor_cauchy() calculates Cauchy correlation with nugget effect", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    corr_0 <- cor_cauchy(x, a = 5, alpha = 0.5, nu = 3)
    corr_0 <- add_nugget(corr_0, nugget = 0.5)
    corr <- cor_cauchy(x, a = 5, alpha = 0.5, nu = 3, nugget = 0.5,
                       is.dist = TRUE)
    expect_equal(corr_0, corr)
})

test_that("cor_cauchy() errors if any(x) < 0", {
    x <- matrix(c(0, -5, 0, -5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(cor_cauchy(x = x, a = 5, alpha = 0.5, is.dist = TRUE))
})
