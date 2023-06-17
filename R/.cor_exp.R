#' Calculate exponential correlation
#'
#' @param c Scale parameter
#' @param gamma Smooth parameter
#' @param x A vector, matrix, or array
#'
#' @keywards internal
#' @export
.cor_exp <- function(c, gamma, x) {
    return(exp(-c * x ^ (2 * gamma)))
}
