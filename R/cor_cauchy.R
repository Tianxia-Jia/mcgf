#' Calculate Cauchy correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param a Smooth parameter, \eqn{a>0}.
#' @param alpha Scale parameter, \eqn{\alpha\in(0, 1]}.
#' @param nu Power parameter, \eqn{\nu>0}. Default is 1.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The Cauchy correlation function with scale parameter \eqn{a} and
#' smooth parameter \eqn{\alpha} has the form
#' \deqn{C(x)=(a|x|^{2\alpha} + 1)^{-\nu}.}
#'
#' @references
#' Gneiting, T., and Schlather, M. (2004). Stochastic Models That Separate
#' Fractal Dimension and the Hurst Effect. SIAM Review, 46(2), 269–282.
.cor_cauchy <- function(x, a, alpha, nu = 1) {
    return((a * abs(x) ^ (2 * alpha) + 1) ^ (-nu))
}

#' Calculate Cauchy correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param a Smooth parameter, \eqn{a>0}.
#' @param alpha Scale parameter, \eqn{\alpha\in(0, 1]}.
#' @param nu Power parameter, \eqn{\nu>0}. Default is 1.
#' @param nugget The nugget effect \eqn{\in[0, 1]}. Used only when `is.dist` is
#' TRUE.
#' @param is.dist Logical; if TRUE, `x` is a distance matrix or an array of
#' distance matrices.
#'
#' @return Correlations with the same dimension as `x`.
#' @export
#'
#' @details
#' The Cauchy correlation function with scale parameter \eqn{a} and
#' smooth parameter \eqn{\alpha} has the form
#' \deqn{C(x)=(1-\text{nugget})(a|x|^{2\alpha} + 1)^{-\nu}+\text{nugget}\cdot
#' \delta_{x=0},} where \eqn{\delta_{x=0}} is 1 when \eqn{x=0} and 0 otherwise.
#'
#' @examples
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_cauchy(x, a = 1, alpha = 0.5)
#'
#' x <- matrix(c(0, 5, 5, 0), nrow = 2)
#' cor_cauchy(x, a = 1, alpha = 0.5, nugget = 0.3, is.dist = TRUE)
#'
#' @references
#' Gneiting, T., and Schlather, M. (2004). Stochastic Models That Separate
#' Fractal Dimension and the Hurst Effect. SIAM Review, 46(2), 269–282.
#'
#' @seealso [cor_exp], [cor_fs], [cor_lagr_tri]
cor_cauchy <- function(x, a, alpha, nu = 1, nugget = 0, is.dist = FALSE) {

    stopifnot(nugget >= 0 && nugget <= 1)
    stopifnot(a > 0)
    stopifnot(alpha > 0 && alpha <= 1)
    stopifnot(nu > 0)

    corr <- .cor_cauchy(a = a, alpha = alpha, nu = nu, x = x)

    if (nugget > 0 & is.dist == F)
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
