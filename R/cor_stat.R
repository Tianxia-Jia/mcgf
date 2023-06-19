#' Calculate general stationary correlation.
#'
#' @param base Matrix or array of correlations.
#' @param dist List of distance matrices or arrays.
#' @param lag_max Maximum lag considered for the Lagrangian term.
#' @param par List of parameters for the Lagrangian term.
#' @param lambda
#' @param asymm_model
#'
#' @return
#' @export
#'
#' @examples
.corr_stat <- function(base, dist, lag_max, par, lambda,
                       asymm_model = c("none", "lagr_tri")){

    asymm_model <- match.arg(asymm_model)

    if (asymm_model == "none") {

        return(base)

    } else if (asymm_model == "lagr_tri") {

        corr_lagr <- cor_lagr_tri(v1 = par$v1, v2 = par$v2, h1 = dist$h1,
                                  h2 = dist$h2, u = lag_max)
        corr <- (1 - lambda) * base + lambda * corr_lagr

        return(corr)
    }
}

corr_stat <- function() {

}
