test_that("cor_fs() calculates fs correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.001, gamma = 0.25)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 5, alpha = 0.1)
    fit_sep <- cor_sep(
        spatial = "exp", temporal = "cauchy",
        par_s = par_s, par_t = par_t, h = h, u = u
    )
    fit_fs <- cor_fs(
        nugget = par_s$nugget, c = par_s$c, gamma = par_s$gamma,
        a = par_t$a, alpha = par_t$alpha, beta = 0, h = h, u = u
    )
    expect_equal(fit_sep, fit_fs)
})

test_that("cor_fs() calculates fs correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.001, gamma = 0.25)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 5, alpha = 0.1)
    fit_sep <- cor_sep(
        spatial = "exp", temporal = "cauchy",
        par_s = par_s, par_t = par_t, h = h, u = u
    )
    fit_fs <- cor_fs(
        nugget = par_s$nugget, c = par_s$c, gamma = par_s$gamma,
        a = par_t$a, alpha = par_t$alpha, beta = 0, h = h, u = u
    )
    expect_equal(fit_sep, fit_fs)
})
