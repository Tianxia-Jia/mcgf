test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 1, alpha = 0.5, nu = 3)
    corr <- cor_sep(
        spatial = "exp", temporal = "cauchy",
        par_s = par_s, par_t = par_t, h = h, u = u
    )

    fit_s <- cor_exp(
        x = h, c = par_s$c, gamma = par_s$gamma,
        nugget = par_s$nugget, is.dist = T
    )
    fit_t <- cor_cauchy(x = u, a = par_t$a, alpha = par_t$alpha, nu = par_t$nu)
    expect_equal(corr, fit_s * fit_t)
})

test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(c = 0.02, gamma = 0.2)
    corr <- cor_sep(
        spatial = "exp", temporal = "exp",
        par_s = par_s, par_t = par_t, h = h, u = u
    )

    fit_s <- cor_exp(
        x = h, c = par_s$c, gamma = par_s$gamma,
        nugget = par_s$nugget, is.dist = T
    )
    fit_t <- cor_exp(x = u, c = par_t$c, gamma = par_t$gamma)
    expect_equal(corr, fit_s * fit_t)
})

test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, a = 1, alpha = 0.5, nu = 5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 5, alpha = 0.8, nu = 3)
    corr <- cor_sep(
        spatial = "cauchy", temporal = "cauchy",
        par_s = par_s, par_t = par_t, h = h, u = u
    )

    fit_s <- cor_cauchy(
        x = h, a = par_s$a, alpha = par_s$alpha, nu = par_s$nu,
        nugget = par_s$nugget, is.dist = T
    )
    fit_t <- cor_cauchy(x = u, a = par_t$a, alpha = par_t$alpha, nu = par_t$nu)
    expect_equal(corr, fit_s * fit_t)
})
