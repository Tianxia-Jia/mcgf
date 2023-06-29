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
#' each matrix slice) of `x` add `nugget`.
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

#' Adjust for nugget effect for correlations
#'
#' @param x A correlation matrix or 3-d array.
#' @param nugget A scalar nugget effect.
#' @param set A correlation matrix or 3-d array of the same dimension as `x`.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' The nugget effect is set to the correlations by first multiplying every
#' entry by \eqn{(1-\text{nugget})}. Then the diagonals (or the diagonals of
#' each matrix slice) of `x` are set to the corresponding diagonals of `set`.
set_nugget <- function(x, nugget, set) {

    if (any(dim(x) != dim(set)))
        stop("dimensions for `x` and `set` must match.")

    corr <- (1 - nugget) * x
    dim_x <- dim(x)

    if (length(dim_x) == 3) {
        for (i in 1:dim_x[3])
            diag(corr[, , i]) <- diag(set[, , i])
    } else if (length(dim_x) == 2) {
        diag(corr) <- diag(set)
    } else {
        stop("invalid dimention for 'x'.")
    }
    return(corr)
}
