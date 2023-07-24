#' Covariance for joint distribution
#'
#' @name cov_joint
#'
#' @param cov Array of covariance matrices.
#'
#' @keywords internal
#' @return A partition of the covariance matrix for the joint distribution of
#' the current value and the past value for a Markov chain Gaussian field.
#'
#' @details
#' The covariance matrix of the joint distribution has the block toeplitz
#' structure. Input `cov` is assumed to be an array of covariance matrices
#' where the \eqn{i}th matrix slice correspond to the \eqn{(i-1)}th time lag. For
#' example, `cov[, , 1]` is the covariance matrix for time lag 0. All matrices
#' in `cov` are used to construct the joint covariance matrix.
#'
#' `cov_joint` gives the joint covariance matrix.
#'
#' `cov_curr` gives the covariance matrix for the current observations.
#'
#' `cov_curr_past` gives partition of the joint covariance for the current and
#' past observations.
#'
#' `cov_past_current` gives partition of the joint covariance for the past and
#' current observations.
#'
#' `cov_past` gives the covariance matrix for the past observations.
cov_joint <- function(cov) {

    n_var <- dim(cov)[1]
    lag_max <- dim(cov)[3] - 1
    n_block <- dim(cov)[3]
    ind_toep <- stats::toeplitz(0:lag_max)

    cov_all <- NULL
    for (i in 1:n_block) {
        ind_u_i <- ind_toep[i, ]
        cov_all <- rbind(cov_all,
                         matrix(cov[, , ind_u_i + 1], nrow = n_var))
    }
    return(cov_all)
}

#' @rdname cov_joint
#' @param horizon Forecast horizon.
cov_curr <- function(cov, horizon) {
    n_var <- dim(cov)[1]
    ind_toep <- stats::toeplitz(0:(horizon - 1))

    cov_mat <- NULL
    for (i in 1:horizon) {
        ind_horizon_i <- ind_toep[i, ]
        cov_mat <- rbind(cov_mat,
                         matrix(cov[, , ind_horizon_i + 1], nrow = n_var))
    }
    return(cov_mat)
}

#' @rdname cov_joint
#' @param lag Time lag.
cov_curr_past <- function(cov, lag) {
    n_var <- dim(cov)[1]
    lag_max <- dim(cov)[3] - 1
    n_block <- dim(cov)[3]
    horizon <- n_block - lag
    ind_toep <- stats::toeplitz(0:lag_max)[1:horizon, -c(1:horizon),
                                           drop = FALSE]

    cov_mat <- NULL
    for (i in 1:horizon) {
        ind_lag_i <- ind_toep[i, ]
        cov_mat <- cbind(cov_mat,
                         t((apply(cov[, , ind_lag_i + 1], 1, c))))
    }
    return(cov_mat)
}

#' @rdname cov_joint
#' @param lag Time lag.
cov_past_curr <- function(cov, lag) {
    n_var <- dim(cov)[1]
    lag_max <- dim(cov)[3] - 1
    n_block <- dim(cov)[3]
    horizon <- n_block - lag
    ind_toep <- stats::toeplitz(0:lag_max)[-c(1:horizon), 1:horizon,
                                           drop = FALSE]

    cov_mat <- NULL
    for (i in 1:horizon) {
        ind_lag_i <- ind_toep[, i]
        cov_mat <- cbind(cov_mat,
                         apply(cov[, , ind_lag_i + 1], 2, c))
    }
    return(cov_mat)
}

#' @rdname cov_joint
#' @param lag Time lag.
cov_past <- function(cov, lag) {
    n_var <- dim(cov)[1]
    ind_toep <- stats::toeplitz(0:(lag - 1))

    cov_mat <- NULL
    for (i in 1:lag) {
        ind_lag_i <- ind_toep[i, ]
        cov_mat <- rbind(cov_mat,
                         matrix(cov[, , ind_lag_i + 1], nrow = n_var))
    }
    return(cov_mat)
}
