#' Simulate regime-switching Markov chain Gaussian field
#'
#' @param N Sample size.
#' @param labels Vector of regime labels of the same length as `N`.
#' @param base_ls List of base model, `sep` or `fs` for now.
#' @param lagrangian_ls List of Lagrangian model, "none" or `lagr_tri` for now.
#' @param par_base_ls List of parameters for the base model.
#' @param par_lagr_ls List of parameters for the Lagrangian model.
#' @param lambda_ls List of weight of the Lagrangian term,
#' \eqn{\lambda\in[0, 1]}.
#' @param dists_ls List of distance matrices or arrays.
#' @param sd_ls List of standard deviation for each location.
#' @param lag_ls List of time lags.
#' @param horizon Forecast horizon, default is 1.
#' @param init Initial samples, default is 0.
#' @param mu_c,mu_p List of means of current and past.
#' @param return_all Logical; if TRUE the joint covariance matrix, arrays of
#' distances and time lag are returned.
#'
#' @keywords internal
#'
#' @return Simulated regime-switching Markov chain Gaussian field with
#' user-specified covariance structure. The simulation is done by kriging.
#' The output data is in space-wide format. Each element in `dists_ls` must
#' contain `h` for symmetric models, and `h1` and `h2` for general stationary
#' models. `horizon` controls forecasting horizon. `init` can be a scalar or a
#' vector of appropriate size. List elements in `sd_ls`, `mu_c_ls`, and
#' `mu_p_ls` must be vectors of appropriate sizes.
.mcgf_rs_sim <- function(N,
                         labels,
                         base_ls,
                         lagrangian_ls,
                         par_base_ls,
                         par_lagr_ls,
                         lambda_ls,
                         dists_ls,
                         sd_ls,
                         lag_ls,
                         horizon = 1,
                         init = 0,
                         mu_c_ls,
                         mu_p_ls,
                         return_all = FALSE) {

    n_regime <- length(unique(labels))

    lag_max_ls <- lapply(lag_ls, function(x) x + horizon - 1)
    n_var <- nrow(dists_ls[[1]]$h)
    n_block_row_ls <- lapply(lag_max_ls, function(x) x + 1)

    u_ls <- lapply(lag_max_ls, function(x) 0:x)
    dim_ar_ls <- lapply(u_ls, function(x) c(n_var, n_var, length(x)))
    h_ar_ls <- Map(function(x, dim) array(x$h, dim = dim),
                   dists_ls, dim_ar_ls)
    u_ar_ls <- Map(function(x, dim)
        array(rep(x, each = n_var * n_var), dim = dim),
        u_ls, dim_ar_ls)

    if (any(lagrangian_ls != "none")) {
        h1_ar_ls <- Map(function(x, dim) array(x["h1"], dim = dim),
                        dists_ls, dim_ar_ls)
        h2_ar_ls <- Map(function(x, dim) array(x["h2"], dim = dim),
                        dists_ls, dim_ar_ls)

        cov_ar_rs <- cor_stat_rs(
            n_regime = n_regime,
            base_ls = base_ls,
            lagrangian_ls = lagrangian_ls,
            par_base_ls = par_base_ls,
            par_lagr_ls = par_lagr_ls,
            lambda_ls = lambda_ls,
            h_ls = h_ls,
            h1_ls = h1_ls,
            h2_ls = h2_ls,
            u_ls = u_ls,
            base_fixed = FALSE
        )

    } else {

        cov_ar_rs <- cor_stat_rs(
            n_regime = n_regime,
            base_ls = base_ls,
            lagrangian_ls = lagrangian_ls,
            par_base_ls = par_base_ls,
            lambda_ls = lambda_ls,
            h_ls = h_ls,
            u_ls = u_ls,
            base_fixed = FALSE
        )
    }

    for (k in 1:n_regime) {
        for (i in 1:dim(cov_ar_rs[[k]])[3]) {
            cov_ar_rs[[k]][, , i] <- cor2cov(cov_ar_rs[[k]][, , i], sd[[k]])
        }
    }




}

mcgf_rs_sim <- function(N,
                        labels,
                        base_ls,
                        lagrangian_ls,
                        par_base_ls,
                        par_lagr_ls,
                        lambda_ls,
                        dists_ls,
                        sd_ls,
                        lag_ls,
                        horizon = 1,
                        init = 0,
                        mu_c_ls = list(0),
                        mu_p_ls = list(0),
                        return_all = FALSE) {

    if (N < horizon)
        stop("'N' must be no less than 'horizon'")

    n_regime <- length(unique(labels))

    if (n_regime == 1)
        cat("Only 1 regime found in 'labels'. Simulating for 1 regime only.")

    for (k in 1:length(dist_ls)) {
        if (is.null(dist_ls[[k]]$h))
            stop("missing 'h' for regime ", k, " in dists_ls.")
    }

    for (k in 1:length(lagrangian_ls)) {

        if (lagrangian_ls[[k]] != "none") {

            if (is.null(dist_ls[[k]]$h1))
                stop("missing 'h1' for regime ", k, " in dists_ls.")
            if (is.null(dist_ls[[k]]$h2))
                stop("missing 'h2' for regime ", k, " in dists_ls.")
        }
    }

    lag_max_ls <- lapply(lag_ls, function(x) x + horizon - 1)
    n_var <- nrow(dists_ls[[1]]$h)
    n_block_row_ls <- lapply(lag_max_ls, function(x) x + 1)

    if (N < horizon + max(unlist(n_block_row_ls))) {
        warning("'N' must be no less than ",
                horizon + max(unlist(n_block_row_ls)))
        N <- horizon + max(unlist(n_block_row_ls))
    }

    if (length(init) == 1) {
        init <- matrix(init, nrow = max(unlist(n_block_row_ls)), ncol = n_var)
    } else {
        if (NROW(init) != max(unlist(n_block_row_ls)) || NCOL(init) !=  n_var)
            stop("dim of 'n_var' must be 1 or ", max(unlist(n_block_row_ls)),
                 " x ", n_var, ".")
    }


    sd_ls <- check_length_ls(sd_ls, n_var, "sd_ls")
    mu_c_ls <- check_length_ls(mu_c_ls, n_var * horizon, "mu_c_ls")
    mu_p_ls <- check_length_ls(mu_p_ls, n_var * horizon, "mu_p_ls")

    res <- .mcgf_rs_sim(
        N = N,
        labels = labels,
        base_ls = base_ls,
        lagrangian_ls = lagrangian_ls,
        par_base_ls = par_base_ls,
        par_lagr_ls = par_lagr_ls,
        lambda_ls = lambda_ls,
        dists_ls = dists_ls,
        sd_ls = sd_ls,
        lag_ls = lag_ls,
        horizon = horizon,
        init = init,
        mu_c_ls = mu_c_ls,
        mu_p_ls = mu_p_ls,
        return_all = return_all
    )
    return(res)
}
