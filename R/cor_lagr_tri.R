#' Calculate Lagrangian correlation of the triangular form
#'
#' @param v1 Prevailing wind, u-component.
#' @param v2 Prevailing wind, v-component.
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
#' \deqn{C(\mathbf{h}, u)=\left(1-\dfrac{1}{2\|\boldsymbol v\|}
#' \left|\dfrac{\mathbf{h}^\top\boldsymbol v}{\|\boldsymbol v\|}-
#' u\|\boldsymbol v\|\right|\right)_+,}
#' where \eqn{x_+=\max(x,0).}
.cor_lagr_tri <- function(v1, v2, h1, h2, u) {
    v_norm <- sqrt(v1 ^ 2 + v2 ^ 2)

    lagr <- 1 - 1 / (2 * v_norm) *
        abs((h1 * v1 + h2 * v2) / v_norm - v_norm * u)
    lagr[lagr < 0] <- 0

    return(lagr)
}

#' Calculate Lagrangian correlation of the triangular form
#'
#' @param v1 Prevailing wind, u-component.
#' @param v2 Prevailing wind, v-component.
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
#' \deqn{C(\mathbf{h}, u)=\left(1-\dfrac{1}{2\|\boldsymbol v\|}
#' \left|\dfrac{\mathbf{h}^\top\boldsymbol v}{\|\boldsymbol v\|}-
#' u\|\boldsymbol v\|\right|\right)_+,}
#' where \eqn{x_+=\max(x,0).}
#'
#' @examples
#' h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
#' h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
#' u <- matrix(0.1, nrow = 2, ncol = 2)
#' cor_lagr_tri(v1 = 5, v = 10, h1 = h1, h2 = h2, u = u)
#'
#' h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
#' h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
#' u <- array(rep(c(0.1, 0.2, 0.3), each = 4), dim = c(2, 2, 3))
#' cor_lagr_tri(v1 = 5, v = 10, h1 = h1, h2 = h2, u = u)
#'
#' @seealso [cor_exp], [cor_cauchy], [cor_fs]
cor_lagr_tri <- function(v1, v2, h1, h2, u) {

    if (!(length(dim(h1)) %in% 2:3))
        stop("'h1' must be a matrix or 3-d array")
    if (any(dim(h1) != dim(h2)))
        stop("'h2' must be of the same dimension as 'h1'")
    if (any(dim(h1) != dim(u)))
        stop("'u' must be of the same dimension as 'h1'")

    corr <- .cor_lagr_tri(v1 = v1, v2 = v2, h1 = h1, h2 = h2, u = u)

    return(corr)
}
