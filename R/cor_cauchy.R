#' Calculate Cauchy correlation
#'
#' @param x A numeric vector, matrix, or array.
#' @param a Smooth parameter, \eqn{a>0}.
#' @param alpha Scale parameter, \eqn{\alpha\in(0, 1]}.
#' @param nu Power parameter, \eqn{\nu>0}. Default is 1.
#' @param nugget The nugget effect \eqn{\in[0, 1]}.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The Cauchy correlation function with scale parameter \eqn{a} and
#' smooth parameter \eqn{\alpha} has the form
#' \deqn{C(x)=(1-\text{nugget})(a|x|^{2\alpha} + 1)^{-\nu}+\text{nugget}\cdot
#' \delta_{x=0},} where \eqn{\delta_{x=0}} is 1 when \eqn{x=0} and 0 otherwise.
#'
#' @references
#' Gneiting, T., and Schlather, M. (2004). Stochastic Models That Separate
#' Fractal Dimension and the Hurst Effect. SIAM Review, 46(2), 269–282.
.cor_cauchy <- function(x, a, alpha, nu = 1, nugget = 0) {
    corr <- (a * abs(x) ^ (2 * alpha) + 1) ^ (-nu)
    if (nugget > 0)
        corr <- add_nugget(x = corr, nugget = nugget)
    return(corr)
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
#' @seealso [cor_exp], [cor_fs], [cor_sep], [cor_lagr_tri], [cor_stat]
cor_cauchy <- function(x, a, alpha, nu = 1, nugget = 0, is.dist = FALSE) {

    if (!is_numeric_scalar(nugget) || nugget < 0 || nugget > 1)
        stop('"nugget" must be in [0, 1].')

    if (!is_numeric_scalar(a) || a <= 0)
        stop('"a" must be positive.')

    if (!is_numeric_scalar(alpha) || alpha <= 0 || alpha > 1)
        stop('"alpha" must be in (0, 1].')

    if (!is_numeric_scalar(nu) || nu <= 0)
        stop('"nu" must be positive.')

    if (nugget > 0 && is.dist == F)
        stop('nugget effect used only when "is.dist = TRUE".')

    if (is.dist) {
        check_dist(x)
        corr <- .cor_cauchy(x = x, a = a, alpha = alpha, nu = nu,
                            nugget = nugget)
    } else {
        corr <- .cor_cauchy(x = x, a = a, alpha = alpha, nu = nu, nugget = 0)
    }
    return(corr)
}
