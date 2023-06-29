#' Check if numeric scalar
#'
#' @param x Input
#'
#' @return Logical.
#' @keywords internal
#'
#' @details
#' Check if input `x` is a numeric  scalar.
is_numeric_scalar <- function(x) {
    return(is.atomic(x) && length(x) == 1L && is.numeric(x) && !is.na(x))
}

#' Check if valid distance
#'
#' @param x Distance matrix or array.
#'
#' @return `x`.
#' @keywords internal
#'
#' @details
#' Check if input `x` is a valid distance vector, matrix or array. It errors if
#' any elements in `x` is negative, or if `x` is not a symmetric matrix or an
#' array of symmetric matrices.
check_dist <- function(x) {
    if (any(is.na(x))) {
        stop('NA or NaN found in "x".')
    }

    if (any(x < 0)) {
        stop('invalid negative distance in "x".')
    }

    if (!is_numeric_scalar(x)) {
        if (is.array(x)) {
            if (!(length(dim(x)) %in% 2:3)) {
                stop('"x" must be a matrix or 3-d array.')
            }

            if (is.matrix(x)) {
                if (!isSymmetric.matrix(x))
                    stop("distance matrix 'x' is not symmetric.")
            } else {
                for (i in 1:dim(x)[3])
                    if (!isSymmetric.matrix(x[, , i]))
                        stop('not all matrix slices in array "x" is symmetric.')
            }
        }
    }
    return(invisible(NULL))
}

#' Check if valid signed distance
#'
#' @param x Distance matrix or array.
#'
#' @return `x`.
#' @keywords internal
#'
#' @details
#' Check if input `x` is a valid signed distance vector, matrix or array. It
#' errors if `x` in absolute value is not a symmetric matrix or an array of
#' symmetric matrices.
check_dist_sign <- function(x) {
    if (any(is.na(x))) {
        stop('NA or NaN found in "x".')
    }

        if (!is_numeric_scalar(x)) {
        if (is.array(x)) {
            if (!(length(dim(x)) %in% 2:3)) {
                stop('"x" must be a matrix or 3-d array.')
            }

            if (is.matrix(x)) {
                if (!isSymmetric.matrix(abs(x)))
                    stop("distance matrix 'x' is not symmetric in absolute
                         values.")
            } else {
                for (i in 1:dim(x)[3])
                    if (!isSymmetric.matrix(abs(x)[, , i]))
                        stop('not all matrix slices in array "x" is symmetric
                             in absolute values.')
            }
        }
    }
    return(invisible(NULL))
}
