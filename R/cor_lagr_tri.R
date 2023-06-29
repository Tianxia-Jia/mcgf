#' Calculate Lagrangian correlation of the triangular form
#'
#' @param v1 Prevailing wind, u-component.
#' @param v2 Prevailing wind, v-component.
#' @param k Scale parameter of \eqn{\|\boldsymbol v\|}, \eqn{k>0}. Default is 2.
#' @param h1 Horizontal distance matrix or array.
#' @param h2 Vertical distance matrix or array, same dimension as `h1`.
#' @param u Time lag, same dimension as `h1`.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `h1`.
#'
#' @details
#' The Lagrangian correlation function of the triangular form with parameters
#' \eqn{\boldsymbol v = (v_1, v_2)^\top} has the form
#' \deqn{C(\mathbf{h}, u)=\left(1-\dfrac{1}{k\|\boldsymbol v\|}
#' \left|\dfrac{\mathbf{h}^\top\boldsymbol v}{\|\boldsymbol v\|}-
#' u\|\boldsymbol v\|\right|\right)_+,}
#' where \eqn{x_+=\max(x,0), \mathbf{h} = (h_1, h_2)^\top}, and \eqn{k > 0} is
#' the scale parameter controls the magnitude of asymmetry in correlation.
.cor_lagr_tri <- function(v1, v2, k = 2, h1, h2, u) {
    v_norm <- sqrt(v1 ^ 2 + v2 ^ 2)

    lagr <- 1 - 1 / (k * v_norm) *
        abs((h1 * v1 + h2 * v2) / v_norm - v_norm * u)
    lagr[lagr < 0] <- 0

    return(lagr)
}

#' Calculate Lagrangian correlation of the triangular form
#'
#' @param v1 Prevailing wind, u-component.
#' @param v2 Prevailing wind, v-component.
#' @param k Scale parameter of \eqn{\|\boldsymbol v\|}, \eqn{k>0}. Default is 2.
#' @param h1 Horizontal distance matrix or array.
#' @param h2 Vertical distance matrix or array, same dimension as `h1`.
#' @param u Time lag, same dimension as `h1`.
#'
#' @return Correlations of the same dimension as `h1`.
#' @export
#'
#' @details
#' The Lagrangian correlation function of the triangular form with parameters
#' \eqn{\boldsymbol v = (v_1, v_2)^\top} has the form
#' \deqn{C(\mathbf{h}, u)=\left(1-\dfrac{1}{k\|\boldsymbol v\|}
#' \left|\dfrac{\mathbf{h}^\top\boldsymbol v}{\|\boldsymbol v\|}-
#' u\|\boldsymbol v\|\right|\right)_+,}
#' where \eqn{x_+=\max(x,0), \mathbf{h} = (h_1, h_2)^\top}, and \eqn{k > 0} is
#' the scale parameter controls the magnitude of asymmetry in correlation.
#'
#' @examples
#' h1 <- matrix(c(0, -5, 5, 0), nrow = 2)
#' h2 <- matrix(c(0, -8, 8, 0), nrow = 2)
#' u <- matrix(0.1, nrow = 2, ncol = 2)
#' cor_lagr_tri(v1 = 5, v2 = 10, h1 = h1, h2 = h2, u = u)
#'
#' h1 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
#' h2 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
#' u <- array(rep(-c(1, 2, 3), each = 4), dim = c(2, 2, 3))
#' cor_lagr_tri(v1 = 10, v2 = 10, h1 = h1, h2 = h2, u = u)
#'
#' @seealso [cor_exp], [cor_cauchy], [cor_fs], [cor_stat]
cor_lagr_tri <- function(v1, v2, k = 2, h1, h2, u) {

    if (!is_numeric_scalar(k) || k <= 0)
        stop('"k" must be positive.')

    if (any(dim(h1) != dim(h2)))
        stop("'h1' must be of the same dimension as 'h2'.")

    if (any(dim(h1) != dim(u)))
        stop("'u' must be of the same dimension as 'h1'.")

    check_dist_sign(h1)
    check_dist_sign(h2)

    corr <- .cor_lagr_tri(v1 = v1, v2 = v2, k = k, h1 = h1, h2 = h2, u = u)

    return(corr)
}
