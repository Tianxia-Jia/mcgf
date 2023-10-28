#' Adjust for nugget effect for correlations
#'
#' @name add_nugget
#'
#' @param x A correlation matrix or 3-d array.
#' @param nugget A scalar nugget effect.
#'
#' @keywords internal
#'
#' @return Correlations of the same dimension as `x`.
#'
#' @details
#' To adjust spatial nugget effect, enery entry of `x` is first multipled by
#' by \eqn{(1-\text{nugget})}; Then `add_nugget` adds `nugget` to the diagonals
#' (or the diagonals of each matrix slice) of `x`, and `set_nugget` set the
#' diagonals (or the diagonals of each matrix slice) to the corresponding
#' diagonals of `set_to`.
add_nugget <- function(x, nugget) {
    corr <- (1 - nugget) * x
    dim_x <- dim(x)

    if (length(dim_x) == 3) {
        for (i in 1:dim_x[3]) {
            diag(corr[, , i]) <- diag(corr[, , i]) + nugget
        }
    } else if (length(dim_x) == 2) {
        diag(corr) <- diag(corr) + nugget
    } else {
        stop("invalid dimention for 'x'.")
    }
    return(corr)
}

#' @rdname add_nugget
#' @param set_to A correlation matrix or 3-d array of the same dimension as `x`.
set_nugget <- function(x, nugget, set_to) {
    if (any(dim(x) != dim(set_to))) {
        stop("dimensions for `x` and `set_to` must match.")
    }

    corr <- (1 - nugget) * x
    dim_x <- dim(x)

    if (length(dim_x) == 3) {
        for (i in 1:dim_x[3]) {
            diag(corr[, , i]) <- diag(set_to[, , i])
        }
    } else if (length(dim_x) == 2) {
        diag(corr) <- diag(set_to)
    } else {
        stop("invalid dimention for `x`.")
    }
    return(corr)
}
