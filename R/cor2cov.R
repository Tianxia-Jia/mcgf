#' Convert correlation to covariance
#'
#' @param V A correlation matrix, usually positive semi-definite.
#' @param sd A vector of standard deviations.
#'
#' @return A correlation matrix.
#' @export
#'
#' @examples
#' V <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' sd <- 1:2
#' cor2cov(V, sd)
cor2cov <- function(V, sd) {

    p <- (d <- dim(V))[1L]

    if (!is.numeric(V) || length(d) != 2L || p != d[2L])
        stop("`V` is not a square numeric matrix", call. = FALSE)

    if (any(V < 0))
        stop("`V` must be non-negative", call. = FALSE)

    stopifnot(dim(V) == c(length(sd), length(sd)))

    sd_mat <- as.matrix(sd)
    return(V * sd_mat %*% t(sd_mat))
}
