#' Calculate general stationary correlation.
#'
#' @param par_base Parameters for the base model (symmetric), used only when
#' `base_fixed = FALSE`.
#' @param par_lagr Parameters for the Lagrangian model.
#' @param base Base model, `sep` or `fs` for now. Or correlation matrix/array.
#' @param lagrangian Lagrangian model, `lagr_tri` for now.
#' @param base_fixed Logical; if TRUE, `base` is the correlation.
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
.corr_stat <- function(par_base,
                       par_lagr,
                       base,
                       lagrangian,
                       base_fixed = FALSE) {
    if (base_fixed) {
        fit_base <- base
    } else {
        fit_base <- do.call(paste0("cor_", base), par_base)
    }

    fit_lagr <- do.call(paste0("cor_", lagrangian))

    stopifnot(dim(fit_base) == dim(fit_lagr))
    corr <- (1 - par_lagr$lambda) * base + par_lagr$lambda * lagrangian

    return(corr)
}

#' Calculate general stationary correlation.
#'
#' @param par_base Parameters for the base model (symmetric), used only when
#' `base_fixed = FALSE`.
#' @param par_lagr Parameters for the Lagrangian model.
#' @param h Euclidean distance matrix or array.
#' @param h1 Horizontal distance matrix or array, same dimension as `h`.
#' @param h2 Vertical distance matrix or array, same dimension as `h`.
#' @param u Time lag, same dimension as `h`.
#' @param base Base model, `sep` or `fs` for now. Or correlation matrix/array.
#' @param lagrangian Lagrangian model, `lagr_tri` for now.
#' @param base_fixed Logical; if TRUE, `base` is the correlation.
#'
#' @return Correlations for the general stationary model. Same dimension of
#' `base` if `base_fixed = FALSE`.
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
corr_stat <- function(par_base,
                      par_lagr,
                      h,
                      h1,
                      h2,
                      u,
                      base = c("sep", "fs"),
                      lagrangian = c("lagr_tri"),
                      base_fixed = FALSE) {
    if (base_fixed) {
        fit_base <- base
    } else {
        fit_base <- do.call(paste0("cor_", base), par_base)
    }

    fit_lagr <- do.call(paste0("cor_", lagrangian))

    stopifnot(dim(fit_base) == dim(fit_lagr))
    corr <- (1 - par_lagr$lambda) * base + par_lagr$lambda * lagrangian

    return(corr)
}
