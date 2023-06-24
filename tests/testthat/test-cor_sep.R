test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 1, alpha = 0.5, nu = 3)
    corr <- cor_sep(par_s = par_s, par_t = par_t, h = h, u = u,
                    spatial = "exp", temporal = "cauchy")

    fit_s <- cor_exp(x = h, c = par_s$c, gamma = par_s$gamma, is.dist = T)
    fit_t <- cor_cauchy(x = u, a = par_t$a, alpha = par_t$alpha, nu = 3)
    fit_sep <- add_nugget(fit_s * fit_t, nugget = par_s$nugget)
    expect_equal(corr, fit_sep)
})

test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(c = 0.02, gamma = 0.2)
    corr <- cor_sep(par_s = par_s, par_t = par_t, h = h, u = u,
                    spatial = "exp", temporal = "exp")

    fit_s <- cor_exp(x = h, c = par_s$c, gamma = par_s$gamma, is.dist = T)
    fit_t <- cor_exp(x = u, c = par_t$c, gamma = par_t$gamma)
    fit_sep <- add_nugget(fit_s * fit_t, nugget = par_s$nugget)
    expect_equal(corr, fit_sep)
})

test_that("cor_sep() calculates separable correlation", {
    h <- array(c(0, 5, 0, 5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, a = 1, alpha = 0.5, nu = 5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 5, alpha = 0.8, nu = 3)
    corr <- cor_sep(par_s = par_s, par_t = par_t, h = h, u = u,
                    spatial = "cauchy", temporal = "cauchy")

    fit_s <- cor_cauchy(x = h, a = par_s$a, alpha = par_s$alpha, nu = par_s$nu,
                        is.dist = T)
    fit_t <- cor_cauchy(x = u, a = par_t$a, alpha = par_t$alpha, nu = par_t$nu)
    fit_sep <- add_nugget(fit_s * fit_t, nugget = par_s$nugget)
    expect_equal(corr, fit_sep)
})

test_that("cor_sep() errors if any(h) < 0", {
    h <- array(c(0, 5, 0, -5, 0, 1, 0, 1, 0), dim = c(3, 3, 3))
    par_s <- list(nugget = 0.5, a = 1, alpha = 0.5, nu = 5)
    u <- array(c(1, 2, 3), dim = c(3, 3, 3))
    par_t <- list(a = 5, alpha = 0.8, nu = 3)
    expect_error(cor_sep(par_s = par_s, par_t = par_t, h = h, u = u,
                         spatial = "cauchy", temporal = "cauchy"))
})
