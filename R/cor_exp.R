#' Calculate exponential correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param c Smooth parameter, \eqn{c>0}.
#' @param gamma Scale parameter, \eqn{\gamma\in(0, 1/2]}. Default is 1/2.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The exponential correlation function with scale parameter \eqn{c} and
#' smooth parameter \eqn{\gamma} has the form
#' \deqn{C(x)=\exp(-c|x|^{2\gamma}).}
#'
#' @references
#' Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
#' Geostatistics. Journal of the Royal Statistical Society. Series C (Applied
#' Statistics), 47(3), 299–350.
.cor_exp <- function(x, c, gamma = 1 / 2) {
    return(exp(-c * abs(x) ^ (2 * gamma)))
}

#' Calculate exponential correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param c Smooth parameter, \eqn{c>0}.
#' @param gamma Scale parameter, \eqn{\gamma\in(0, 1/2]}. Default is 1/2.
#' @param nugget The nugget effect \eqn{\in[0, 1]}. Used only when `is.dist` is
#' TRUE.
#' @param is.dist Logical; if TRUE, `x` is a distance matrix or an array of
#' distance matrices.
#'
#' @return Correlations with the same dimension as `x`.
#' @export
#'
#' @details
#' The exponential correlation function with scale parameter \eqn{c}
#' and smooth parameter \eqn{\gamma} has the form
#' \deqn{C(x)=(1-\text{nugget})\exp(-c|x|^{2\gamma})+\text{nugget}\cdot
#' \delta_{x=0},} where \eqn{\delta_{x=0}} is 1 when \eqn{x=0} and 0 otherwise.
#'
#' @examples
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_exp(x, c = 0.01, gamma = 0.5)
#'
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_exp(x, c = 0.01, gamma = 0.5, nugget = 0.3, is.dist = TRUE)
#'
#' @references
#' Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
#' Geostatistics. Journal of the Royal Statistical Society. Series C (Applied
#' Statistics), 47(3), 299–350.
#'
#' @seealso [cor_cauchy], [cor_fs], [cor_sep], [cor_lagr_tri]
cor_exp <- function(x, c, gamma = 1 / 2, nugget = 0, is.dist = FALSE) {

    stopifnot(nugget >= 0 & nugget <= 1)
    stopifnot(c > 0)
    stopifnot(gamma > 0 & gamma <= 1 / 2)

    corr <- .cor_exp(c = c, gamma = gamma, x = x)

    if (nugget > 0 && is.dist == F)
        stop("nugget effect used only when 'is.dist = TRUE'.")

    if (is.dist) {
        if (any(x < 0))
            stop("invalid negative distance in 'x'.")
        else if (!(length(dim(x)) %in% 2:3))
            stop("'x' must be a matrix or 3-d array")

        if (nugget > 0)
            corr <- add_nugget(x = corr, nugget = nugget)
    }
    return(corr)
}
