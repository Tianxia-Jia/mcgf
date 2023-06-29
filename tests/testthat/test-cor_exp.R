test_that("cor_exp() calculates exponential correlation", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    corr <- cor_exp(x, c = 0.01, gamma = 0.5)
    expect_equal(corr, .cor_exp(x = x, c = 0.01, gamma = 0.5))
})

test_that("cor_exp() errors if nugget > 0 but is.dist = FALSE", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(cor_exp(x, c = 0.01, gamma = 0.5, nugget = 0.5,
                         is.dist = FALSE))
})

test_that("cor_cauchy() errors if invalid distance", {
    x <- matrix(c(-1, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    expect_error(cor_cauchy(x = x, a = 5, alpha = 0.5, nu = 3, is.dist = TRUE))
})


test_that("cor_exp() calculates exponential correlation with nugget effect", {
    x <- matrix(c(0, 5, 0, 5, 0, 1, 0, 1, 0), nrow = 3)
    corr_0 <- cor_exp(x, c = 0.01, gamma = 0.5)
    corr_0 <- add_nugget(corr_0, nugget = 0.5)
    corr <- cor_exp(x, c = 0.01, gamma = 0.5, nugget = 0.5, is.dist = TRUE)
    expect_equal(corr_0, corr)
})
