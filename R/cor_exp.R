#' Calculate exponential correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param c Smooth parameter, \eqn{c>0}.
#' @param gamma Scale parameter, \eqn{\gamma\in(0, 1/2]}. Default is 1/2.
#' @param nugget The nugget effect \eqn{\in[0, 1]}.
#'
#' @keywords internal
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The exponential correlation function with scale parameter \eqn{c}
#' and smooth parameter \eqn{\gamma} has the form
#' \deqn{C(x)=(1-\text{nugget})\exp(-c|x|^{2\gamma})+\text{nugget}\cdot
#' \delta_{x=0},} where \eqn{\delta_{x=0}} is 1 when \eqn{x=0} and 0 otherwise.
#'
#' @references
#' Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
#' Geostatistics. Journal of the Royal Statistical Society. Series C (Applied
#' Statistics), 47(3), 299–350.
.cor_exp <- function(x, c, gamma = 1 / 2, nugget = 0) {
    corr <- exp(-c * abs(x)^(2 * gamma))
    if (nugget > 0) {
        corr <- add_nugget(x = corr, nugget = nugget)
    }
    return(corr)
}

#' Calculate exponential correlation
#'
#' @inherit .cor_exp return params details references
#'
#' @param is.dist Logical; if TRUE, `x` is a distance matrix or an array of
#' distance matrices.
#'
#' @export
#' @examples
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_exp(x, c = 0.01, gamma = 0.5)
#'
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_exp(x, c = 0.01, gamma = 0.5, nugget = 0.3, is.dist = TRUE)
#'
#' @family {correlation functions}
cor_exp <- function(x, c, gamma = 1 / 2, nugget = 0, is.dist = FALSE) {
    if (!is_numeric_scalar(nugget) || nugget < 0 || nugget > 1) {
        stop("`nugget` must be in [0, 1].", call. = FALSE)
    }

    if (!is_numeric_scalar(c) || c <= 0) {
        stop("`c` must be positive.", call. = FALSE)
    }

    if (!is_numeric_scalar(gamma) || gamma <= 0 || gamma > 1 / 2) {
        stop("`gamma` must be in (0, 1/2].", call. = FALSE)
    }

    corr <- .cor_exp(c = c, gamma = gamma, x = x)

    if (nugget > 0 && is.dist == F) {
        stop("nugget effect used only when `is.dist = TRUE`.", call. = FALSE)
    }

    if (is.dist) {
        check_dist(x = x)
        corr <- .cor_exp(x = x, c = c, gamma = gamma, nugget = nugget)
    } else {
        corr <- .cor_exp(x = x, c = c, gamma = gamma, nugget = 0)
    }
    return(corr)
}
