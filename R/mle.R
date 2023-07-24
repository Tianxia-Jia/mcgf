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

    cov_mat_curr <- cov_curr(cov = fitted, horizon = horizon)
    cov_mat_curr_past <- cov_curr_past(cov = fitted, lag = lag)
    cov_mat_past_curr <- cov_past_curr(cov = fitted, lag = lag)
    cov_mat_past <- cov_past(cov = fitted, lag = lag)

    cov_mat_past_inv <- tryCatch({
        solve(cov_mat_past)
    },
    error = function(e) {
        MASS::ginv(cov_mat_past)
    })

    weights <- cov_mat_curr_past %*% cov_mat_past_inv
    Sigma_c <- cov_mat_curr - weights %*% t(cov_mat_curr_past)

    det_Sigma_c <- det(Sigma_c)

    if (is.na(det_Sigma_c) | det_Sigma_c < 0) {
        return(Inf)
    } else {
        x_ts <- stats::embed(as.matrix(x), n_lag)

        mu_c <- t(tcrossprod(weights, x_ts[, -c(1:(n_var * horizon))]))
        mu_diff <- x_ts[, 1:(n_var * horizon)] - mu_c

        Sigma_c_inv <- tryCatch({
            solve(Sigma_c)
        },
        error = function(e) {
            MASS::ginv(Sigma_c)
        })

        llike <- - nrow(x_ts) * log(det_Sigma_c) -
            sum(apply(mu_diff, 1, function(x, y) t(x) %*% y %*% x, Sigma_c_inv))

        return(-llike)
    }
}
