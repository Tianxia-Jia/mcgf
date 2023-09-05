#' Compute the objective for wls method
#'
#' @param par Parameters of `cor_fn`.
#' @param cor_fn Correlation function.
#' @param cor_emp Empirical correlations.
#' @param par_fixed Fixed parameters of `cor_fn`.
#'
#' @keywords internal
#' @return The objective of weighted least squares.
obj_wls <- function(par, cor_fn, cor_emp, par_fixed) {

    fitted <- do.call(cor_fn, c(par, par_fixed))
    summand <- ((cor_emp - fitted) / (1 - fitted)) ^ 2
    summand[is.infinite(summand)] <- NA
    wls <- sum(summand, na.rm = T)
    return(wls)
}

#' Title
#'
#' @param par Parameters of `cor_fn`.
#' @param cor_fn Correlation function
#' @param x An `mcgf` or `mcgf_rs` object
#' @param lag Time lag.
#' @param par_fixed Fixed parameters of `cor_fn`.
#'
#' @keywords internal
#' @return The objective of maximum likelihood: the additive inverse of
#' log-likelihood.
obj_mle <- function(par, cor_fn, x, lag, par_fixed) {

    sds <- sds(x)
    n_var <- length(sds)

    fitted <- do.call(cor_fn, c(par, par_fixed))
    n_lag <- dim(fitted)[3]

    for (i in 1:n_lag) {
        fitted[, , i] <- .cor2cov(V = fitted[, , i], sd = sds)
    }

    lag_max <- n_lag - 1
    horizon <- n_lag - lag

    new_cov_par <- cov_par(cov = fitted, horizon = horizon)
    det_cov_curr <- det(new_cov_par$cov_curr)

    if (is.na(det_cov_curr) | det_cov_curr < 0) {
        return(0)
    } else {

        mle_label <- attr(x, "mle_label", exact = TRUE)
        x_ts <- stats::embed(as.matrix(x), n_lag)

        if (!is.null(mle_label) && is.mcgf_rs(x)) {
            label <- attr(x, "label", exact = TRUE)
            ind_n <- label == mle_label
            ind_n <- ind_n[n_lag:nrow(x)]
            x_ts <- x_ts[ind_n, ]
        }

        mu_c <- t(tcrossprod(new_cov_par$weights,
                             x_ts[, -c(1:(n_var * horizon))]))
        mu_diff <- x_ts[, 1:(n_var * horizon)] - mu_c

        cov_curr_inv <- mat_inv(new_cov_par$cov_curr)

        llike <- -nrow(x_ts) * log(det_cov_curr) -
            sum(apply(mu_diff, 1, function(x, y)
                t(x) %*% y %*% x, cov_curr_inv))

        if (is.infinite(llike)) {
            return(0)
        } else {
            return(-llike)
        }
    }
}

#' Optimization for wls method
#'
#' @param par_init Initial values for parameters to be optimized.
#' @param method Parameter estimation method. "wls" or "mle".
#' @param optim_fn Optimization function.
#' @param cor_fn Correlation function.
#' @param par_fixed Fixed parameters of `cor_fn`.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param ... Additional arguments passed to `optim_fn`, `obj_wls` or `obj_mle`.
#'
#' @keywords internal
#' @return A list outputted from optimization functions of `optim_fn`.
estimate <- function(par_init, method, optim_fn, cor_fn, par_fixed, lower,
                     upper, ...) {

    obj_fn <- switch(method,
                     wls = obj_wls,
                     mle = obj_mle)

    args <- list(par_init,
                 obj_fn,
                 cor_fn = cor_fn,
                 par_fixed = par_fixed,
                 lower = lower,
                 upper = upper,
                 ...)

    if (optim_fn == "optim") args <- c(args, method = "L-BFGS-B")

    return(do.call(optim_fn, args))
}
