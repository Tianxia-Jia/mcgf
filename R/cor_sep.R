#' Calculate correlation for separable model
#'
#' @param par_s Parameters for the pure spatial model. Nugget effect supported.
#' @param par_t Parameters for the pure temporal model.
#' @param spatial Pure spatial model, `exp` or `cauchy` for now.
#' @param temporal Pure temporal model, `exp` or `cauchy` for now.
#'
#' @keywords internal
#' @return Correlations for separable model.
#'
#' @details
#' The separable model is the product of a pure spatial model and pure temporal
#' model. Now only exponential and Cauchy correlation models are available.
.cor_sep <- function(par_s, par_t, spatial, temporal) {

    fit_s <- do.call(paste0("cor_", spatial), par_s)
    fit_t <- do.call(paste0("cor_", temporal), par_t)
    fit_sep <- fit_s * fit_t

    return(fit_sep)
}

#' Calculate correlation for separable model
#'
#' @param par_s Parameters for the pure spatial model. Nugget effect supported.
#' @param par_t Parameters for the pure temporal model.
#' @param h Euclidean distance matrix or array.
#' @param u Time lag, same dimension as `h`.
#' @param spatial Pure spatial model, `exp` or `cauchy` for now.
#' @param temporal Pure temporal model, `exp` or `cauchy` for now.
#'
#' @return Correlations of the same dimension as u
#' @export
#'
#' @details
#' The separable model is the product of a pure spatial model and pure temporal
#' model. Now only exponential and Cauchy correlation models are available.
#'
#' @examples
#' h <- matrix(c(0, 5, 5, 0), nrow = 2)
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' u <- matrix(0, nrow = 2, ncol = 2)
#' par_t <- list(a = 1, alpha = 0.5)
#' cor_sep(par_s = par_s, par_t = par_t, h = h, u = u, spatial = "exp",
#'         temporal = "cauchy")
#'
#' h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
#' par_t <- list(a = 1, alpha = 0.5)
#' cor_sep(par_s = par_s, par_t = par_t, h = h, u = u, spatial = "exp",
#'         temporal = "cauchy")
#'
#' @references
#' Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
#' Space–Time Data, Journal of the American Statistical Association, 97:458,
#' 590-600.
#'
#' @seealso [cor_exp], [cor_cauchy], [cor_fs], [cor_lagr_tri], [cor_stat]
cor_sep <- function(par_s,
                    par_t,
                    h,
                    u,
                    spatial = c("exp", "cauchy"),
                    temporal = c("exp", "cauchy")) {

    spatial <- match.arg(spatial)
    temporal <- match.arg(temporal)

    if (any(h < 0))
        stop("invalid negative distance in 'h'.")

    if (!(length(dim(h)) %in% 2:3))
        stop("'h' must be a matrix or 3-d array")

    if (length(dim(h)) == 2 && !isSymmetric.matrix(h))
        stop("distance matrix 'h' is not symmetric.")

    if (length(dim(h)) == 3) {
        for (i in 1:dim(h)[3])
            if (!isSymmetric.matrix(h[,,i]))
                stop("distance array 'h' is not symmetric")
    }

    if (any(dim(h) != dim(u)))
        stop("'u' must be of the same dimension as 'h' in 'dist'")

    nugget <- par_s[["nugget"]]
    par_s[["nugget"]] <- NULL
    par_s <- append(par_s, list(x = h, is.dist = TRUE))
    par_t <- append(par_t, list(x = u, is.dist = FALSE))

    corr <- .cor_sep(par_s = par_s, par_t = par_t, spatial = spatial,
                     temporal = temporal)
    if (nugget > 0)
        corr <- add_nugget(x = corr, nugget = nugget)
    return(corr)
}
