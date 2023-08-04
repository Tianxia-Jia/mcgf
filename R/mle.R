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
        return(Inf)
    } else {
        x_ts <- stats::embed(as.matrix(x), n_lag)

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
