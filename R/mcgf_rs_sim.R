#' Simulate regime-switching Markov chain Gaussian field
#'
#' @param N Sample size.
#' @param labels Vector regime labels of the same length as `N`.
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
#' models. `horizon` controls forecasting horizon. `sd`, `init`, `mu_c`, and
#' `mu_p` must be vectors of appropriate sizes.
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



    u <- 0:lag_max
    dim_ar <- c(n_var, n_var, length(u))
    h_ar <- array(dists$h, dim = dim_ar)
    u_ar <- array(rep(u, each = n_var * n_var), dim = dim_ar)



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

        if (lagrangian_ls[[k]] == "lagr_tri") {

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


    for (k in length(sd_ls)) {
        if (length(sd_ls[[k]]) == 1) {
            sd_ls[[k]] <- rep(sd_ls[[k]], n_var)
        } else {
            if (length(sd_ls[[k]]) != n_var)
                stop("length of 'sd_ls' must be 1 or ", n_var, " for regime ", k)
        }
    }

    if (length(init) == 1) {
        init <- matrix(init, nrow = max(unlist(n_block_row_ls)), ncol = n_var)
    } else {
        if (NROW(init) != max(unlist(n_block_row_ls)) || NCOL(init) !=  n_var)
            stop("dim of 'n_var' must be 1 or ", max(unlist(n_block_row_ls)),
                 " x ", n_var)
    }

    for (k in length(mu_c_ls)) {
        if (length(mu_c_ls[[k]]) == 1) {
            mu_c_ls[[k]] <- rep(mu_c_ls[[k]], n_var * horizon)
        } else {
            if (length(mu_c_ls[[k]]) != n_var * horizon)
                stop("length of 'mu_c_ls' must be 1 or ", n_var * horizon,
                     " for regime ", k)
        }
    }

    for (k in length(mu_p_ls)) {
        if (length(mu_p_ls[[k]]) == 1) {
            mu_p_ls[[k]] <- rep(mu_p_ls[[k]], n_var * lag_ls[[k]])
        } else {
            if (length(mu_p_ls[[k]]) != n_var * lag_ls[[k]])
                stop("length of 'mu_p_ls' must be 1 or ", n_var * lag_ls[[k]],
                     " for regime ", k)
        }
    }

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
