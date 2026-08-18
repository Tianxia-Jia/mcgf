make_test_data <- function(n = 60L) {
    t <- seq_len(n)
    data.frame(
        S1 = sin(t / 3) + 0.01 * t,
        S2 = cos(t / 4) + sin(t / 7),
        S3 = sin(t / 5) - cos(t / 6)
    )
}

make_test_locations <- function() {
    matrix(
        c(0, 0, 2, 0, 1, 1.5),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(c("S1", "S2", "S3"), c("x", "y"))
    )
}

make_test_dists <- function() {
    find_dists(make_test_locations(), longlat = FALSE)
}

make_test_mcgf <- function(n = 60L) {
    suppressMessages(mcgf(
        make_test_data(n),
        locations = make_test_locations(),
        longlat = FALSE
    ))
}

make_test_mcgf_rs <- function(n = 60L) {
    suppressMessages(
        mcgf_rs(
            make_test_data(n),
            locations = make_test_locations(),
            label = rep(c("A", "B"), length.out = n),
            longlat = FALSE
        )
    )
}

make_mock_base_fit <- function(x, lag = 2L, horizon = 1L) {
    par <- c(
        c = 0.2,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.5
    )
    par_fixed <- c(nugget = 0)

    lag_max <- lag + horizon - 1L
    hu <- mcgf:::to_ar(dists(x)$h, lag_max = lag_max)
    cor_base <- mcgf:::..cor_sep(
        nugget = par_fixed[["nugget"]],
        c = par[["c"]],
        gamma = par[["gamma"]],
        a = par[["a"]],
        alpha = par[["alpha"]],
        h = hu$h_ar,
        u = hu$u_ar
    )

    list(
        model = "sep",
        method = "wls",
        optim_fn = "nlminb",
        lag = lag,
        horizon = horizon,
        scale_time = 1,
        fit = list(
            par = unname(par),
            objective = 0,
            convergence = 0
        ),
        par_names = names(par),
        par_fixed = par_fixed,
        dists_base = NULL,
        dots = list(),
        par_base = as.list(c(par, par_fixed)),
        cor_base = cor_base,
        rs = FALSE
    )
}

make_mock_lagr_fit <- function(x) {
    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    scale_time <- attr(x, "scale_time", exact = TRUE)
    lag_max <- lag + horizon - 1L

    d <- dists(x)
    h1 <- mcgf:::to_ar(d$h1, lag_max = lag_max, u = FALSE)
    h2 <- mcgf:::to_ar(d$h2, lag_max = lag_max, u = FALSE)
    u <- mcgf:::to_ar(d$h, lag_max = lag_max)$u_ar / scale_time

    par <- c(
        lambda = 0.15,
        v1 = 1.5,
        v2 = 0.5
    )
    par_fixed <- c(k = 2)

    cor_lagr <- mcgf:::..cor_stat(
        cor_base = attr(x, "base_res", exact = TRUE)$cor_base,
        lagrangian = "lagr_tri",
        lambda = par[["lambda"]],
        v1 = par[["v1"]],
        v2 = par[["v2"]],
        k = par_fixed[["k"]],
        h1 = h1,
        h2 = h2,
        u = u
    )

    list(
        model = "lagr_tri",
        method = "wls",
        optim_fn = "nlminb",
        fit = list(
            par = unname(par),
            objective = 0,
            convergence = 0
        ),
        par_names = names(par),
        par_fixed = par_fixed,
        dists_base = FALSE,
        dists_lagr = NULL,
        dots = list(),
        par_lagr = as.list(c(par, par_fixed)),
        cor_lagr = cor_lagr,
        rs = FALSE
    )
}
