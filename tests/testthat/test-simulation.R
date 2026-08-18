test_that("mcgf_sim() produces the expected dimensions for a base-only model", {
    d <- make_test_dists()
    par_base <- list(
        par_s = list(
            nugget = 0.05,
            c = 0.1,
            gamma = 0.5
        ),
        par_t = list(a = 0.2, alpha = 0.5)
    )

    set.seed(1)
    x <- mcgf_sim(
        N = 20,
        base = "sep",
        lagrangian = "none",
        par_base = par_base,
        par_lagr = NULL,
        lambda = 0,
        dists = d,
        lag = 1
    )

    expect_equal(ncol(x), 3)
    expect_equal(nrow(x), 22)
    expect_equal(colnames(x), colnames(d$h))
    expect_true(all(is.finite(x)))
})

test_that("mcgf_sim() can return the covariance construction inputs", {
    d <- make_test_dists()
    par_base <- list(
        par_s = list(
            nugget = 0.05,
            c = 0.1,
            gamma = 0.5
        ),
        par_t = list(a = 0.2, alpha = 0.5)
    )

    set.seed(2)
    out <- mcgf_sim(
        N = 10,
        base = "sep",
        lagrangian = "none",
        par_base = par_base,
        par_lagr = NULL,
        lambda = 0,
        dists = d,
        lag = 1,
        return_all = TRUE
    )

    expect_named(out, c("X", "par"))
    expect_true(is.matrix(out$par$cov_mat))
    expect_equal(dim(out$par$u), c(3, 3, 2))
})

test_that("mcgf_rs_sim() produces regime labels and data columns", {
    d <- make_test_dists()
    par_base <- list(
        par_s = list(
            nugget = 0.05,
            c = 0.1,
            gamma = 0.5
        ),
        par_t = list(a = 0.2, alpha = 0.5)
    )
    label <- rep(c(1, 2), length.out = 20)

    set.seed(3)
    x <- mcgf_rs_sim(
        N = 20,
        label = label,
        base_ls = list("sep"),
        lagrangian_ls = list("none"),
        par_base_ls = list(par_base),
        par_lagr_ls = list(NULL),
        lambda_ls = list(0),
        dists_ls = list(d),
        lag_ls = list(1, 1)
    )

    expect_equal(ncol(x), 4)
    expect_equal(colnames(x), c("regime", colnames(d$h)))
    expect_equal(unname(tail(x[, "regime"], length(label))), label)
})
