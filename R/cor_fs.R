#' Calculate correlation for fully symmetric model
#'
#' @param nugget The nugget effect \eqn{\in[0, 1]}.
#' @param c Scale parameter of `cor_exp`, \eqn{c>0}.
#' @param gamma Smooth parameter of `cor_exp`, \eqn{\gamma\in(0, 1/2]}.
#' @param a Scale parameter of `cor_cauchy`, \eqn{a>0}.
#' @param alpha Smooth parameter of `cor_cauchy`, \eqn{\alpha\in(0, 1]}.
#' @param beta Interaction parameter, \eqn{\beta\in[0, 1]}.
#' @param h Euclidean distance matrix or array.
#' @param u Time lag, same dimension as `h`.
#'
#' @keywords internal
#' @return Correlations of the same dimension as `h` and `u`.
#'
#' @details
#' The fully symmetric correlation function with interaction parameter
#' \eqn{\beta} has the form
#' \deqn{C(\mathbf{h}, u)=\dfrac{1}{(a|u|^{2\alpha} + 1)}
#' \left((1-\text{nugget})\exp\left(\dfrac{-c\|\mathbf{h}\|^{2\gamma}}
#' {(a|u|^{2\alpha}+1)^{\beta\gamma}}\right)+
#' \text{nugget}\cdot \delta_{\mathbf{h}=\boldsymbol 0}\right),}
#' where \eqn{\|\cdot\|} is the Euclidean distance, and \eqn{\delta_{x=0}} is 1
#' when \eqn{x=0} and 0 otherwise. Here \eqn{\mathbf{h}\in\mathbb{R}^2} and
#' \eqn{u\in\mathbb{R}}. By default `beta = 0` and it reduces to the separable
#' model.
#'
#' @references
#' Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
#' Space–Time Data, Journal of the American Statistical Association, 97:458,
#' 590-600.
.cor_fs <- function(nugget, c, gamma = 1 / 2, a, alpha, beta = 0, h, u) {
    c_cauchy <- .cor_cauchy(a = a, alpha = alpha, nu = 1, x = u)
    c_exp <- .cor_exp(c = c, gamma = gamma, x = h)
    corr <- c_cauchy * c_exp ^ (c_cauchy ^ (beta * gamma))
    corr <- set_nugget(x = corr, nugget = nugget, set = c_cauchy)
    return(corr)
}

#' Calculate correlation for fully symmetric model
#'
#' @inherit .cor_fs params details return references
#'
#' @export
#' @examples
#' h <- matrix(c(0, 5, 5, 0), nrow = 2)
#' u <- matrix(0, nrow = 2, ncol = 2)
#' cor_fs(c = 0.01, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#'
#' h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
#' u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
#' cor_fs(c = 0.01, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#'
#' @family {correlation functions}
cor_fs <- function(nugget = 0, c, gamma = 1 / 2, a, alpha, beta = 0, h, u) {

    if (!is_numeric_scalar(nugget) || nugget < 0 || nugget > 1)
        stop('"nugget" must be in [0, 1].')

    if (!is_numeric_scalar(c) || c <= 0)
        stop('"c" must be positive.')

    if (!is_numeric_scalar(gamma) || gamma <= 0 || gamma > 1 / 2)
        stop('"gamma" must be in (0, 1/2].')

    if (!is_numeric_scalar(a) || a <= 0)
        stop('"a" must be positive.')

    if (!is_numeric_scalar(alpha) || alpha <= 0 || alpha > 1)
        stop('"alpha" must be in (0, 1].')

    check_dist(x = h)

    if (any(dim(h) != dim(u)))
        stop("'u' must be of the same dimension as 'h'.")

    corr <- .cor_fs(nugget = nugget, c = c, gamma = gamma, a = a, alpha = alpha,
                    beta = beta, h = h, u = u)
    return(corr)
}
