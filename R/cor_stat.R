#' Calculate general stationary correlation.
#'
#' @param par_base Parameters for the base model (symmetric), used only when
#' `base_fixed = FALSE`.
#' @param par_lagr Parameters for the Lagrangian model.
#' @param base Base model, `sep` or `fs` for now. Or correlation matrix/array.
#' @param lagrangian Lagrangian model, `lagr_tri` for now.
#' @param base_fixed Logical; if TRUE, `base` is the correlation.
#' @param lambda Weight of the Lagrangian term, \eqn{\lambda\in[0, 1]}.
#'
#' @keywords internal
#' @return Correlations for the general stationary model. Same dimension of
#' `base` if `base_fixed = FALSE`.
#'
#' @details The general station model, a convex combination of a base model
#' and a Lagrangian model, has the form
#' \deqn{C(\mathbf{h}, u)=(1-\lambda)C_{\text{Base}}(\mathbf{h}, u)+
#' \lambda C_{\text{Lagr}}(\mathbf{h}, u),}
#' where \eqn{\lambda} is the weight of the Lagrangian term.
#'
#' If `base_fixed = TRUE`, then the correlation is of the form
#' \deqn{C(\mathbf{h}, u)=(1-\lambda)\text{Base}+
#' \lambda C_{\text{Lagr}}(\mathbf{h}, u).}
.corr_stat <- function(par_base = NULL,
                       par_lagr,
                       base,
                       lagrangian,
                       lambda = lambda,
                       base_fixed = FALSE) {
    fit_lagr <- do.call(paste0("cor_", lagrangian), par_lagr)

    if (base_fixed) {
        fit_base <- base
        stopifnot(dim(fit_base) == dim(fit_lagr))
    } else {
        fit_base <- do.call(paste0("cor_", base), par_base)
    }
    corr <- (1 - lambda) * fit_base + lambda * fit_lagr

    return(corr)
}

#' Calculate general stationary correlation.
#'
#' @param par_base Parameters for the base model (symmetric), used only when
#' `base_fixed = FALSE`.
#' @param par_lagr Parameters for the Lagrangian model.
#' @param h Euclidean distance matrix or array, used only when
#' `base_fixed = FALSE`.
#' @param h1 Horizontal distance matrix or array, same dimension as `h`.
#' @param h2 Vertical distance matrix or array, same dimension as `h`.
#' @param u Time lag, same dimension as `h`.
#' @param lambda Weight of the Lagrangian term, \eqn{\lambda\in[0, 1]}.
#' @param base Base model, `sep` or `fs` for now. Or correlation matrix/array.
#' @param lagrangian Lagrangian model, `lagr_tri` for now.
#' @param base_fixed Logical; if TRUE, `base` is the correlation.
#'
#' @return Correlations for the general stationary model. Same dimension of
#' `base` if `base_fixed = TRUE`.
#' @export
#'
#' @details The general station model, a convex combination of a base model
#' and a Lagrangian model, has the form
#' \deqn{C(\mathbf{h}, u)=(1-\lambda)C_{\text{Base}}(\mathbf{h}, u)+
#' \lambda C_{\text{Lagr}}(\mathbf{h}, u),}
#' where \eqn{\lambda} is the weight of the Lagrangian term.
#'
#' If `base_fixed = TRUE`, then the correlation is of the form
#' \deqn{C(\mathbf{h}, u)=(1-\lambda)\text{Base}+
#' \lambda C_{\text{Lagr}}(\mathbf{h}, u).}
#'
#' @examples
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' par_t <- list(a = 1, alpha = 0.5)
#' par_base <- list(par_s = par_s, par_t = par_t)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
#' h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
#' h <- sqrt(h1 ^ 2 + h2 ^ 2)
#' u <- matrix(0.1, nrow = 2, ncol = 2)
#' corr_stat(par_base = par_base, par_lagr = par_lagr, h = h, h1 = h1, h2 = h2,
#'           u = u, lambda = 0.8, base = "sep", lagrangian = "lagr_tri")
#'
#' h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
#' h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
#' h <- sqrt(h1 ^ 2 + h2 ^ 2)
#' u <- array(rep(c(0.1, 0.2, 0.3), each = 4), dim = c(2, 2, 3))
#' fit_base <- cor_fs(nugget = 0.5, c = 0.01, gamma = 0.5, a = 1, alpha = 0.5,
#'                    beta = 0.0, h = h, u = u)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' corr_stat(par_lagr = par_lagr, h1 = h1, h2 = h2, u = u, lambda = 0.8,
#'           base = fit_base, lagrangian = "lagr_tri", base_fixed = TRUE)
#'
#' @seealso [cor_exp], [cor_cauchy], [cor_sep], [cor_lagr_tri], [cor_stat]
corr_stat <- function(par_base,
                      par_lagr,
                      h,
                      h1,
                      h2,
                      u,
                      lambda,
                      base = c("sep", "fs"),
                      lagrangian = c("lagr_tri"),
                      base_fixed = FALSE) {

    lagrangian <- match.arg(lagrangian)
    par_lagr <- append(par_lagr,
                       list(h1 = h1, h2 = h2, u = u))

    if (base_fixed) {
        corr <- .corr_stat(
            base = base,
            par_lagr = par_lagr,
            lagrangian = lagrangian,
            lambda = lambda,
            base_fixed = TRUE
        )

    } else {
        base <- match.arg(base)
        if (base == "sep") {
            par_base <- list(
                par_s = par_base$par_s,
                par_t = par_base$par_t,
                h = h,
                u = u,
                spatial = "exp",
                temporal = "cauchy"
            )
        } else if (base == "fs") {
            par_base <- append(par_base, list(h = h, u = u))
        }
        corr <- .corr_stat(
            par_base = par_base,
            par_lagr = par_lagr,
            base = base,
            lagrangian = lagrangian,
            lambda = lambda,
            base_fixed = FALSE
        )
    }

    return(corr)
}
