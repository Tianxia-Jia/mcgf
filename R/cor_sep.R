#' Calculate correlation for separable model
#'
#' @param spatial Pure spatial model, `exp` or `cauchy` for now.
#' @param temporal Pure temporal model, `exp` or `cauchy` for now.
#' @param par_s Parameters for the pure spatial model. Nugget effect supported.
#' @param par_t Parameters for the pure temporal model.
#'
#' @keywords internal
#' @return Correlations for separable model.
#'
#' @details
#' The separable model is the product of a pure temporal model, \eqn{C_T(u)},
#' and a pure spatial model, \eqn{C_S(\mathbf{h})}. It is of the form
#' \deqn{C(\mathbf{h}, u)=C_{T}(u)
#' \left[(1-\text{nugget})C_{S}(\mathbf{h})+\text{nugget}
#' \delta_{\mathbf{h}=0}\right],}
#' where \eqn{\delta_{x=0}} is 1 when \eqn{x=0} and 0 otherwise. Here
#' \eqn{\mathbf{h}\in\mathbb{R}^2} and \eqn{u\in\mathbb{R}}. Now only
#' exponential and Cauchy correlation models are available.
#'
#' @references
#' Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
#' Space–Time Data, Journal of the American Statistical Association, 97:458,
#' 590-600.
.cor_sep <- function(spatial, temporal, par_s, par_t) {

    fit_s <- do.call(paste0("cor_", spatial), par_s)
    fit_t <- do.call(paste0("cor_", temporal), par_t)
    fit_sep <- fit_s * fit_t

    return(fit_sep)
}

#' Calculate correlation for separable model
#'
#' @inherit .cor_sep params details references
#'
#' @param h Euclidean distance matrix or array.
#' @param u Time lag, same dimension as `h`.
#'
#' @return Correlations of the same dimension as `h` and `u`.
#' @export
#'
#' @examples
#' h <- matrix(c(0, 5, 5, 0), nrow = 2)
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' u <- matrix(0, nrow = 2, ncol = 2)
#' par_t <- list(a = 1, alpha = 0.5)
#' cor_sep(spatial = "exp", temporal = "cauchy", par_s = par_s, par_t = par_t,
#'         h = h, u = u)
#'
#' h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
#' par_t <- list(a = 1, alpha = 0.5)
#' cor_sep(spatial = "exp", temporal = "cauchy", par_s = par_s, par_t = par_t,
#'         h = h, u = u)
#'
#' @family {correlation functions}
cor_sep <- function(spatial = c("exp", "cauchy"),
                    temporal = c("exp", "cauchy"),
                    par_s,
                    par_t,
                    h,
                    u) {

    spatial <- match.arg(spatial)
    temporal <- match.arg(temporal)

    if (is.null(dim(h)) != is.null(dim(u)) || any(dim(h) != dim(u)))
        stop("'u' must be of the same dimension as 'h'.", call. = FALSE)

    par_s <- append(par_s, list(x = h, is.dist = TRUE))
    par_t <- append(par_t, list(x = u, is.dist = FALSE))

    corr <- .cor_sep(spatial = spatial, temporal = temporal, par_s = par_s,
                     par_t = par_t)
    return(corr)
}
