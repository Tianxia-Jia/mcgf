#' Calculate Lagrangian correlation of the exponential form
#'
#' @param v1 Prevailing wind, u-component.
#' @param v2 Prevailing wind, v-component.
#' @param k Scale parameter of \eqn{\|\boldsymbol v\|}, \eqn{k>0}. Default is 2.
#' @param h1 Horizontal distance matrix or array.
#' @param h2 Vertical distance matrix or array, same dimension as `h1`.
#' @param u Time lag, same dimension as `h1`.
#'
#' @keywords internal
#' @return Correlations of the same dimension as `h1`.
#'
#' @details
#' The Lagrangian correlation function of the exponential form with parameters
#' \eqn{\boldsymbol v = (v_1, v_2)^\top\in\mathbb{R}^2} has the form
#' \deqn{C(\mathbf{h}, u)=\exp\left(-\dfrac{1}{k\|\boldsymbol v\|}
#' \left\|\mathbf{h}-u\boldsymbol v\right\|\right),}
#' where \eqn{\|\cdot\|} is the Euclidean distance,
#' \eqn{\mathbf{h} = (\mathrm{h}_1, \mathrm{h}_2)^\top\in\mathbb{R}^2},
#' and \eqn{k > 0} is the scale parameter controlling the magnitude of
#' asymmetry in correlation.
#'
#' @references
#' Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
#' Geostatistics. Journal of the Royal Statistical Society. Series C (Applied
#' Statistics), 47(3), 299–350.
.cor_lagr_exp <- function(v1, v2, k = 2, h1, h2, u) {
    v_norm <- sqrt(v1^2 + v2^2)

    h_vu_norm <- sqrt((h1 - v1 * u)^2 + (h2 - v2 * u)^2)

    lagr <- exp(-1 / (k * v_norm) * h_vu_norm)
    return(lagr)
}

#' Calculate Lagrangian correlation of the exponential form
#'
#' @inherit .cor_lagr_exp params details return references
#'
#' @export
#' @examples
#' h1 <- matrix(c(0, -5, 5, 0), nrow = 2)
#' h2 <- matrix(c(0, -8, 8, 0), nrow = 2)
#' u <- matrix(0.1, nrow = 2, ncol = 2)
#' cor_lagr_exp(v1 = 5, v2 = 10, h1 = h1, h2 = h2, u = u)
#'
#' h1 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
#' h2 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
#' u <- array(rep(-c(1, 2, 3), each = 4), dim = c(2, 2, 3))
#' cor_lagr_exp(v1 = 10, v2 = 10, h1 = h1, h2 = h2, u = u)
#'
#' @family correlation functions
cor_lagr_exp <- function(v1, v2, k = 2, h1, h2, u) {
    if (!is_numeric_scalar(k) || k <= 0) {
        stop("`k` must be positive.", call. = FALSE)
    }

    if (any(dim(h1) != dim(h2))) {
        stop("`h2` must be of the same dimension as `h1`.", call. = FALSE)
    }

    if (any(dim(h1) != dim(u))) {
        stop("`u` must be of the same dimension as `h1`.", call. = FALSE)
    }

    check_dist_sign(h1, name = "h1")
    check_dist_sign(h2, name = "h2")

    corr <- .cor_lagr_exp(v1 = v1, v2 = v2, k = k, h1 = h1, h2 = h2, u = u)

    return(corr)
}
