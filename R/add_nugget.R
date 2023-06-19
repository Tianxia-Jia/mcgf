#' Adjust for nugget effect for correlations
#'
#' @param x A correlation matrix or 3-d array.
#' @param nugget A scalar nugget effect.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The nugget effect is added to the correlations by first multiplying every
#' entry by \eqn{(1-\text{nugget})}. Then the diagonals (or the diagonals of
#' each matrix slice) of `x` are set to 1 based on the assumption that the
#' correlation is 1 when `x` is 0.
#'
add_nugget <- function(x, nugget) {
    corr <- (1 - nugget) * x
    dim_x <- dim(x)

    if (length(dim_x) == 3) {
        for (i in 1:dim_x[3])
            diag(corr[, , i]) <- diag(corr[, , i]) + nugget
    } else if (length(dim_x) == 2) {
        diag(corr) <- diag(corr) + nugget
    } else {
        stop("invalid dimention for 'x'.")
    }
    return(corr)
}
