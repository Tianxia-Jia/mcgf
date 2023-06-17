#' Calculate cauchy correlation
#'
#' @param a Scale parameter
#' @param alpha Smooth parameter
#' @param x A vector, matrix, or array
#'
#' @keywards internal
#' @export
.cor_cauchy <- function(a, alpha, x) {
    return((a * abs(x) ^ (2 * alpha) + 1) ^ (-1))
}
