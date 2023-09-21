#' Covariance for joint distribution
#'
#' @name cov_joint
#'
#' @param cov Array of covariance matrices.
#'
#' @keywords internal
#' @return The joint covariance matrix for the joint distribution of the current
#' values and the past values for a Markov chain Gaussian field.
#'
#' @details
#' The covariance matrix of the joint distribution has the block toeplitz
#' structure. Input `cov` is assumed to be an array of cross-covariance matrices
#' where the \eqn{i}th matrix slice correspond to the \eqn{(i-1)}th time lag.
#' For example, `cov[, , 1]` is the cross-covariance matrix for time lag 0. All
#' matrices in `cov` are used to construct the joint covariance matrix.
#'
#' `cov_par` gives weights and covariance matrix for the current values..
cov_joint <- function(cov) {
    cov_t <- aperm(cov, c(2, 1, 3))
    n_var <- dim(cov)[1]
    lag_max <- dim(cov)[3] - 1
    n_block <- dim(cov)[3]
    ind_toep <- stats::toeplitz(0:lag_max)

    cov_all <- matrix(cov[, , ind_toep[1, ] + 1], nrow = n_var)

    for (i in 2:n_block) {
        ind_u_i <- ind_toep[i, ][-c(1:(i - 1))]
        ind_u_t_i <- ind_toep[i, ][1:(i - 1)]
        cov_i <- cbind(
            matrix(cov_t[, , ind_u_t_i + 1], nrow = n_var),
            matrix(cov[, , ind_u_i + 1], nrow = n_var)
        )

        cov_all <- rbind(cov_all, cov_i)
    }
    return(cov_all)
}

#' Find inverse of a symmetric positive definite matrix
#'
#' @param x A symmetric and positive definite matrix.
#'
#' @keywords internal
#' @return Inverse of x.
mat_inv <- function(x) {
    x_inv <- tryCatch(
        {
            chol2inv(chol(x))
        },
        error = function(e) {
            tryCatch(
                {
                    solve(x)
                },
                error = function(e) {
                    MASS::ginv(x)
                }
            )
        }
    )
    return(x_inv)
}

#' @rdname cov_joint
#'
#' @param horizon Forecast horizon, default is 1.
#' @param n_var Number of locations.
#' @param joint Logical; True if `cov` is the joint covariance matrix.
#'
#' @keywords internal
cov_par <- function(cov, horizon = 1, n_var, joint = FALSE) {
    if (joint) {
        if (missing(n_var)) {
            stop("missing `n_var`.", call. = FALSE)
        }
        cov_mat_joint <- cov
    } else {
        cov_mat_joint <- cov_joint(cov = cov)
        n_var <- dim(cov)[1]
    }

    ind_curr <- 1:(n_var * horizon)

    cov_mat_curr <- cov_mat_joint[ind_curr, ind_curr]
    cov_mat_curr_past <- cov_mat_joint[ind_curr, -ind_curr]
    cov_mat_past <- cov_mat_joint[-ind_curr, -ind_curr]
    cov_mat_past_curr <- cov_mat_joint[-ind_curr, ind_curr]

    cov_mat_past_inv <- mat_inv(cov_mat_past)

    weights <- cov_mat_curr_past %*% cov_mat_past_inv
    cov_curr <- cov_mat_curr - weights %*% cov_mat_past_curr

    return(list(weights = weights, cov_curr = cov_curr))
}
