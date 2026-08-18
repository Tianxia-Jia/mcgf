test_that("cor_sep() is the product of its spatial and temporal components", {
    h <- array(c(0, 1, 1, 0), dim = c(2, 2, 3))
    u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))

    par_s <- list(
        nugget = 0.1,
        c = 0.2,
        gamma = 0.5
    )
    par_t <- list(a = 0.4, alpha = 0.5)

    got <- cor_sep(
        spatial = "exp",
        temporal = "cauchy",
        par_s = par_s,
        par_t = par_t,
        h = h,
        u = u
    )

    expected <- cor_exp(
        h,
        c = par_s$c,
        gamma = par_s$gamma,
        nugget = par_s$nugget,
        is.dist = TRUE
    ) * cor_cauchy(u, a = par_t$a, alpha = par_t$alpha)

    expect_equal(got, expected)
})

test_that("cor_fs() reduces to the separable model when beta = 0", {
    h <- array(c(0, 1, 1, 0), dim = c(2, 2, 3))
    u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))

    fs <- cor_fs(
        nugget = 0.1,
        c = 0.2,
        gamma = 0.5,
        a = 0.4,
        alpha = 0.5,
        beta = 0,
        h = h,
        u = u
    )

    sep <- cor_sep(
        spatial = "exp",
        temporal = "cauchy",
        par_s = list(
            nugget = 0.1,
            c = 0.2,
            gamma = 0.5
        ),
        par_t = list(a = 0.4, alpha = 0.5),
        h = h,
        u = u
    )

    expect_equal(fs, sep)
})

test_that("cor_stat() returns the base model when lagrangian = none", {
    h <- array(c(0, 1, 1, 0), dim = c(2, 2, 2))
    u <- array(rep(0:1, each = 4), dim = c(2, 2, 2))
    par_base <- list(
        nugget = 0,
        c = 0.2,
        gamma = 0.5,
        a = 0.4,
        alpha = 0.5,
        beta = 0.2
    )

    expect_equal(
        cor_stat(
            base = "fs",
            lagrangian = "none",
            par_base = par_base,
            h = h,
            u = u
        ),
        do.call(cor_fs, c(par_base, list(
            h = h, u = u
        )))
    )
})

test_that("cor_stat() forms the requested convex combination", {
    h1 <- array(c(0, -1, 1, 0), dim = c(2, 2, 2))
    h2 <- array(c(0, -2, 2, 0), dim = c(2, 2, 2))
    h <- sqrt(h1^2 + h2^2)
    u <- array(rep(0:1, each = 4), dim = c(2, 2, 2))

    par_base <- list(
        nugget = 0,
        c = 0.2,
        gamma = 0.5,
        a = 0.4,
        alpha = 0.5,
        beta = 0
    )
    par_lagr <- list(v1 = 2, v2 = 1, k = 3)

    base <- do.call(cor_fs, c(par_base, list(h = h, u = u)))
    lagr <- do.call(cor_lagr_exp, c(par_lagr, list(
        h1 = h1, h2 = h2, u = u
    )))

    expect_equal(
        cor_stat(
            base = "fs",
            lagrangian = "lagr_exp",
            par_base = par_base,
            par_lagr = par_lagr,
            lambda = 0.25,
            h = h,
            h1 = h1,
            h2 = h2,
            u = u
        ),
        0.75 * base + 0.25 * lagr
    )
})

test_that("cor_stat_rs() agrees with regime-wise cor_stat() calls", {
    h1 <- array(c(0, -1, 1, 0), dim = c(2, 2, 2))
    h2 <- array(c(0, -2, 2, 0), dim = c(2, 2, 2))
    h <- sqrt(h1^2 + h2^2)
    u <- array(rep(0:1, each = 4), dim = c(2, 2, 2))

    par_base <- list(
        nugget = 0,
        c = 0.2,
        gamma = 0.5,
        a = 0.4,
        alpha = 0.5,
        beta = 0
    )
    par_lagr <- list(v1 = 2, v2 = 1, k = 3)

    got <- cor_stat_rs(
        n_regime = 2,
        base_ls = list("fs"),
        lagrangian_ls = list("none", "lagr_tri"),
        par_base_ls = list(par_base),
        par_lagr_ls = list(NULL, par_lagr),
        lambda_ls = list(0, 0.2),
        h_ls = list(h),
        h1_ls = list(NULL, h1),
        h2_ls = list(NULL, h2),
        u_ls = list(u)
    )

    expected1 <- cor_stat(
        base = "fs",
        lagrangian = "none",
        par_base = par_base,
        h = h,
        u = u
    )
    expected2 <- cor_stat(
        base = "fs",
        lagrangian = "lagr_tri",
        par_base = par_base,
        par_lagr = par_lagr,
        lambda = 0.2,
        h = h,
        h1 = h1,
        h2 = h2,
        u = u
    )

    expect_equal(got[[1]], expected1)
    expect_equal(got[[2]], expected2)
})
