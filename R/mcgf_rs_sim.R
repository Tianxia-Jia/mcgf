#' Simulate regime-switching Markov chain Gaussian field
#'
#' @param N Sample size.
#' @param label Vector of regime labels of the same length as `N`.
#' @param base_ls List of base model, `sep` or `fs` for now.
#' @param lagrangian_ls List of Lagrangian model, "none" or `lagr_tri` for now.
#' @param par_base_ls List of parameters for the base model.
#' @param par_lagr_ls List of parameters for the Lagrangian model.
#' @param lambda_ls List of weight of the Lagrangian term,
#' \eqn{\lambda\in[0, 1]}.
#' @param dists_ls List of distance matrices or arrays.
#' @param sd_ls List of standard deviation for each location.
#' @param lag_ls List of time lags.
#' @param scale_time Scale of time unit, default is 1. Elements in `lag_ls` are
#' divided by `scale_time`.
#' @param init Initial samples, default is 0.
#' @param mu_c_ls,mu_p_ls List of means of current and past.
#' @param return_all Logical; if TRUE the joint covariance matrix, arrays of
#' distances and time lag are returned.
#'
#' @keywords internal
#'
#' @return Simulated regime-switching Markov chain Gaussian field with
#' user-specified covariance structures. The simulation is done by kriging.
#' The output data is in space-wide format. Each element in `dists_ls` must
#' contain `h` for symmetric models, and `h1` and `h2` for general stationary
#' models. `init` can be a scalar or a vector of appropriate size.
#' List elements in `sd_ls`, `mu_c_ls`, and `mu_p_ls` must be vectors of
#' appropriate sizes.
.mcgf_rs_sim <- function(N,
                         label,
                         base_ls,
                         lagrangian_ls,
                         par_base_ls,
                         par_lagr_ls,
                         lambda_ls,
                         dists_ls,
                         sd_ls,
                         lag_ls,
                         scale_time = 1,
                         init = 0,
                         mu_c_ls,
                         mu_p_ls,
                         return_all = FALSE) {
    n_regime <- length(unique(label))
    regime <- sort(unique(label))

    lag_max_ls <- lag_ls
    n_var <- nrow(dists_ls[[1]]$h)
    n_block_row_ls <- lapply(lag_max_ls, function(x) x + 1)

    u_ls <- lapply(lag_max_ls, function(x) (0:x) / scale_time)
    dim_ar_ls <- lapply(u_ls, function(x) c(n_var, n_var, length(x)))
    h_ar_ls <- Map(
        function(x, dim) array(x$h, dim = dim),
        dists_ls, dim_ar_ls
    )
    u_ar_ls <- Map(
        function(x, dim) {
            array(rep(x, each = n_var * n_var), dim = dim)
        },
        u_ls, dim_ar_ls
    )

    if (any(lagrangian_ls != "none")) {
        h1_ar_ls <- Map(
            function(x, dim) array(x[["h1"]], dim = dim),
            dists_ls, dim_ar_ls
        )
        h2_ar_ls <- Map(
            function(x, dim) array(x[["h2"]], dim = dim),
            dists_ls, dim_ar_ls
        )

        cov_ar_rs <- cor_stat_rs(
            n_regime = n_regime,
            base_ls = base_ls,
            lagrangian_ls = lagrangian_ls,
            par_base_ls = par_base_ls,
            par_lagr_ls = par_lagr_ls,
            lambda_ls = lambda_ls,
            h_ls = h_ar_ls,
            h1_ls = h1_ar_ls,
            h2_ls = h2_ar_ls,
            u_ls = u_ar_ls,
            base_fixed = FALSE
        )
    } else {
        cov_ar_rs <- cor_stat_rs(
            n_regime = n_regime,
            base_ls = base_ls,
            lagrangian_ls = lagrangian_ls,
            par_base_ls = par_base_ls,
            h_ls = h_ar_ls,
            u_ls = u_ar_ls,
            base_fixed = FALSE
        )
    }

    for (k in 1:n_regime) {
        for (i in 1:dim(cov_ar_rs[[k]])[3]) {
            cov_ar_rs[[k]][, , i] <- cor2cov(cov_ar_rs[[k]][, , i], sd_ls[[k]])
        }
    }

    X_cov_par <- lapply(cov_ar_rs, cov_par)
    new_cov_chol <- lapply(X_cov_par, function(x) chol(x$cov_curr))

    X <- init
    for (n in 1:N) {
        regime_n <- which(label[n] == regime)

        X_past <- stats::embed(
            utils::tail(X, lag_max_ls[[regime_n]]),
            lag_max_ls[[regime_n]]
        )
        X_new_mean <- mu_c_ls[[regime_n]] +
            tcrossprod(
                X_cov_par[[regime_n]]$weights,
                X_past - mu_p_ls[[regime_n]]
            )

        # X_new <- mvnfast::rmvn(1, X_new_mean, X_cov_par[[regime_n]]$cov_curr)
        X_new <- crossprod(new_cov_chol[[regime_n]], stats::rnorm(length(X_new_mean)))
        X_new <- matrix(X_new + X_new_mean, ncol = n_var, byrow = T)
        X <- rbind(X, X_new)
    }

    rownames(X) <- 1:nrow(X)
    colnames(X) <- colnames(dists_ls[[1]]$h)
    X <- cbind(regime = c(rep(NA, NROW(init)), label), X)

    if (return_all) {
        cov_mat_joint_ls <- lapply(cov_ar_rs, cov_joint)
        par <- list(
            cov_mat_ls = cov_mat_joint_ls,
            dists_ls = list(h = h_ar_ls),
            u = u_ar_ls
        )

        if (any(lagrangian_ls != "none")) {
            par$dists_ls <- list(h = h_ar_ls, h1 = h1_ar_ls, h2 = h2_ar_ls)
        }
        return(list(X = X, par = par))
    } else {
        return(X = X)
    }
}

#' @inherit .mcgf_rs_sim title params details return
#'
#' @export
#' @examples
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' par_t <- list(a = 1, alpha = 0.5)
#' par_base <- list(par_s = par_s, par_t = par_t)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
#' h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
#' h <- sqrt(h1^2 + h2^2)
#' dists <- list(h = h, h1 = h1, h2 = h2)
#'
#' set.seed(123)
#' label <- sample(1:2, 1000, replace = TRUE)
#' X <- mcgf_rs_sim(
#'     N = 1000,
#'     label = label,
#'     base_ls = list("sep"),
#'     lagrangian_ls = list("none", "lagr_tri"),
#'     lambda_ls = list(0, 0.5),
#'     par_base_ls = list(par_base),
#'     par_lagr_ls = list(NULL, par_lagr),
#'     dists_ls = list(dists, dists)
#' )
#' # plot.ts(X[, -1])
#'
#' @family simulations of Markov chain Gaussian fields
mcgf_rs_sim <- function(N,
                        label,
                        base_ls,
                        lagrangian_ls,
                        par_base_ls,
                        par_lagr_ls,
                        lambda_ls,
                        dists_ls,
                        sd_ls,
                        lag_ls,
                        scale_time = 1,
                        init = 0,
                        mu_c_ls = list(0),
                        mu_p_ls = list(0),
                        return_all = FALSE) {
    n_regime <- length(unique(label))

    if (n_regime == 1) {
        message(
            "Only 1 regime found in `label`. ",
            "Simulating for 1 regime only."
        )
    }

    for (k in 1:length(dists_ls)) {
        if (is.null(dists_ls[[k]]$h)) {
            stop("missing 'h' for regime ", k, " in `dists_ls`.", call. = FALSE)
        }
    }

    for (k in 1:length(lagrangian_ls)) {
        if (lagrangian_ls[[k]] != "none") {
            if (is.null(dists_ls[[k]]$h1)) {
                stop("missing 'h1' for regime ", k, " in `dists_ls`.",
                    call. = FALSE
                )
            }
            if (is.null(dists_ls[[k]]$h2)) {
                stop("missing 'h2' for regime ", k, " in `dists_ls`.",
                    call. = FALSE
                )
            }
        }
    }

    if (missing(sd_ls)) {
        sd_ls <- rep(list(1), n_regime)
    }

    if (missing(lag_ls)) {
        lag_ls <- rep(list(1), n_regime)
    }

    lag_max_ls <- lag_ls
    n_var <- nrow(dists_ls[[1]]$h)
    n_block_row_ls <- lapply(lag_max_ls, function(x) x + 1)

    if (N < 1 + max(unlist(n_block_row_ls))) {
        warning(
            "'N' must be no less than ",
            1 + max(unlist(n_block_row_ls))
        )
        N <- 1 + max(unlist(n_block_row_ls))
    }

    if (length(init) == 1) {
        init <- matrix(init, nrow = max(unlist(n_block_row_ls)), ncol = n_var)
    } else {
        if (NROW(init) != max(unlist(n_block_row_ls)) || NCOL(init) != n_var) {
            stop("dim of `init` must be 1 or ", max(unlist(n_block_row_ls)),
                " x ", n_var, ".",
                call. = FALSE
            )
        }
    }

    sd_ls <- check_length_ls(sd_ls, n_var, "sd_ls")
    mu_c_ls <- check_length_ls(mu_c_ls, n_var, "mu_c_ls")
    mu_p_ls <- check_length_ls(mu_p_ls, n_var, "mu_p_ls")

    if (length(sd_ls) != n_regime) {
        sd_ls <- rep(sd_ls, n_regime)
    }
    if (length(mu_c_ls) != n_regime) {
        mu_c_ls <- rep(mu_c_ls, n_regime)
    }
    if (length(mu_p_ls) != n_regime) {
        mu_p_ls <- rep(mu_p_ls, n_regime)
    }

    res <- .mcgf_rs_sim(
        N = N,
        label = label,
        base_ls = base_ls,
        lagrangian_ls = lagrangian_ls,
        par_base_ls = par_base_ls,
        par_lagr_ls = par_lagr_ls,
        lambda_ls = lambda_ls,
        dists_ls = dists_ls,
        sd_ls = sd_ls,
        lag_ls = lag_ls,
        scale_time = scale_time,
        init = init,
        mu_c_ls = mu_c_ls,
        mu_p_ls = mu_p_ls,
        return_all = return_all
    )
    return(res)
}
