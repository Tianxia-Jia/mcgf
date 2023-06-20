#' Calculate correlation for fully symmetric model
#'
#' @param c Scale parameter of `cor_exp`, \eqn{c>0}.
#' @param gamma Smooth parameter of `cor_exp`, \eqn{\gamma\in(0, 1/2]}.
#' @param a Scale parameter of `cor_cauchy`, \eqn{a>0}.
#' @param alpha Smooth parameter of `cor_cauchy`, \eqn{\alpha\in(0, 1]}.
#' @param beta Interaction parameter, \eqn{\beta\in[0, 1]}.
#' @param h Euclidean distance matrix or array.
#' @param u Time lag, same dimension as `h`.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The fully symmetric correlation function with interaction parameter
#' \eqn{\beta} has the form
#' \deqn{\dfrac{1}{(a|u|^{2\alpha} + 1)}
#' \exp\left(\dfrac{-c\|h\|^{2\gamma}}
#' {(a|u|^{2\alpha}+1)^{\beta\gamma}}\right),}
#'
#' where \eqn{\|\cdot\|} is the Euclidean distance. By default `beta = 0` and
#' it reduce to the separable model.
#'
#' @references
#' Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
#' Space–Time Data, Journal of the American Statistical Association, 97:458,
#' 590-600.
.cor_fs <- function(c, gamma = 1 / 2, a, alpha, beta = 0, h, u) {
    c_cauchy <- .cor_cauchy(a = a, alpha = alpha, nu = 1, x = u)
    c_exp <- .cor_exp(c = c, gamma = gamma, x = h)

    corr <- c_cauchy * c_exp ^ (c_cauchy ^ (beta * gamma))

    return(corr)
}

#' Calculate correlation for fully symmetric model
#'
#' @param c Scale parameter of `cor_exp`, \eqn{c>0}.
#' @param gamma Smooth parameter of `cor_exp`, \eqn{\gamma\in(0, 1/2]}.
#' @param a Scale parameter of `cor_cauchy`, \eqn{a>0}.
#' @param alpha Smooth parameter of `cor_cauchy`, \eqn{\alpha\in(0, 1]}.
#' @param beta Interaction parameter, \eqn{\beta\in[0, 1]}.
#' @param h Euclidean distance matrix or array.
#' @param u Time lag, same dimension as `h`.
#' @param nugget The nugget effect \eqn{\in[0, 1]}.
#'
#' @return Correlations of the same dimension as `h` and `u`.
#' @export
#'
#' @details
#' The fully symmetric correlation function with interaction parameter
#' \eqn{\beta} has the form
#' \deqn{C(\mathbf{h}, u)=\dfrac{\text{nugget}}{(a|u|^{2\alpha} + 1)}
#' \exp\left(\dfrac{-c\|\mathbf{h}\|^{2\gamma}}
#' {(a|u|^{2\alpha}+1)^{\beta\gamma}}\right)+
#' \text{nugget}\cdot\delta_{\|\mathbf{h}\|=0},}
#'
#' where \eqn{\|\cdot\|} is the Euclidean distance, and \eqn{\delta_{x=0}} is 1
#' when \eqn{x=0} and 0 otherwise. Here \eqn{\mathbf{h}\in\mathbb{R}^2} and
#' \eqn{u\in\mathbb{R}}. By default `beta = 0` and it reduce to the separable
#' model.
#'
#' @examples
#' h <- matrix(c(0, 5, 5, 0), nrow = 2)
#' u <- matrix(0, nrow = 2, ncol = 2)
#' cor_fs(c = 1, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#'
#' h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
#' u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
#' cor_fs(c = 1, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#'
#' @references
#' Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
#' Space–Time Data, Journal of the American Statistical Association, 97:458,
#' 590-600.
#'
#' @seealso [cor_exp], [cor_cauchy], [cor_sep], [cor_lagr_tri]
cor_fs <- function(nugget = 0, c, gamma = 1 / 2, a, alpha, beta = 0, h, u) {

    stopifnot(nugget >= 0 & nugget <= 1)
    stopifnot(c > 0)
    stopifnot(gamma > 0 & gamma <= 1 / 2)
    stopifnot(a > 0)
    stopifnot(alpha > 0 & alpha <= 1)

    if (any(h < 0))
        stop("invalid negative distance in 'h'.")
    if (!(length(dim(h)) %in% 2:3))
        stop("'h' must be a matrix or 3-d array")

    if (any(dim(h) != dim(u)))
        stop("'u' must be of the same dimension as 'h'")

    corr <- .cor_fs(c = c, gamma = gamma, a = a, alpha = alpha, beta = beta,
                    h = h, u = u)
    if (nugget > 0)
        corr <- add_nugget(x = corr, nugget = nugget)

    return(corr)
}
