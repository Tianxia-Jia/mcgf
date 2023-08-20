.cor2cov <- function(V, sd) {

    sd_mat <- as.matrix(sd)
    return(V * sd_mat %*% t(sd_mat))
}

#' Convert correlation to covariance
#'
#' @name cor2cov
#'
#' @param V A correlation matrix, usually positive semi-definite.
#' @param sd A vector of standard deviations.
#' @param empirical Logical; TRUE if V is empirical correlation.
#'
#' @return A correlation matrix.
#' @export
#'
#' @details
#' `cor2cov` converts a matrix. `cor2cov_ar` converts an 3-D array.
#'
#'
#' @examples
#' V <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' sd <- 1:2
#' cor2cov(V, sd)
#'
#' V_ar <- array(c(1, 0.5, 0.5, 1), dim = c(2, 2, 2))
#' cor2cov_ar(V_ar, sd)
cor2cov <- function(V, sd, empirical = FALSE) {

    p <- (d <- dim(V))[1L]

    if (!is.numeric(V) || length(d) != 2L || p != d[2L])
        stop("`V` is not a square numeric matrix", call. = FALSE)

    if (!empirical && any(V < 0))
        stop("`V` must be non-negative", call. = FALSE)

    stopifnot(dim(V) == c(length(sd), length(sd)))

    sd_mat <- as.matrix(sd)
    return(.cor2cov(V = V, sd = sd))
}

#' @rdname cor2cov
#' @export
cor2cov_ar <- function(V, sd, empirical = FALSE) {

    for (i in 1:dim(V)[3]) {
        V[, , i] <- cor2cov(V[, , i], sd, empirical = empirical)
    }
    return(V)
}
