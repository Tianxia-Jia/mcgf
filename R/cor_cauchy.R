#' Calculate Cauchy correlation
#'
#' @param nugget The nugget effect in \eqn{\in[0, 1]}. Effective only when is.dist is TRUE.
#' @param a Smooth parameter, a>0
#' @param alpha Scale parameter, \eqn{\alpha\in(0, 1]}.
#' @param x A vector, matrix, or array.
#' @param is.dist Logical; if TRUE, x is a distance matrix or an array of distance matrices.
#'
#' @return sf
#' @export
#'
#' @examples
#' x <- matrix(c(0, 5, 5, 0))
#' cor_cauchu(x, a = 1, alpha = 0.5)
#'
#'
cor_cauchy <- function(nugget = 0, a, alpha, x, is.dist = F) {

    stopifnot(nugget >=0 && nugget <= 1)
    stopifnot(a > 0)
    stopifnot(alpha > 0 && alpha <= 1)

    corr <- .cor_cauchy(a = a, alpha = alpha, x = x)

    if (is.dist == T & nugget > 0) corr <- add_nugget(nugget, corr)

    return(corr)
}
