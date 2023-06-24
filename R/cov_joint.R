#' Covariance for joint distribution
#'
#' @param cov Array of covariance matrices.
#'
#' @return A covariance matrix for the joint distribution of the current value
#' and the past value for a Markov chain Gaussian field.
#' @export
#'
#' @details
#' The covariance matrix of the joint distribution has the block Toeplitz
#' structure. Input `cov` is assumed to be an array of covariance matrices
#' where the ith matrix slice correspond to the i-1 time lag. For example,
#' `cov[, , 1]` is the covariance matrix for time lag 0. All matrices in `cov`
#' are used to construct the joint covariance matrix.
#'
#' @examples
#' 1
cov_joint <- function(cov) {

    n_var <- dim(cov)[1]
    lag_max <- dim(cov)[3] - 1
    n_block <- dim(cov)[3]
    ind_toep <- toeplitz(0:lag_max)

    cov_all <- NULL
    for (i in 1:n_block) {
        ind_u_i <- ind_toep[i, ]
        cov_all <- rbind(cov_all,
                         matrix(cov[, , ind_u_i + 1], nrow = n_var))
    }
    return(cov_all)
}
