#' Check if numeric scalar
#'
#' @param x Input
#'
#' @return Logical.
#' @keywords internal
#'
#' @details
#' Check if `x` is a numeric scalar.
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
#' Check if `x` is a valid distance vector, matrix or array. It errors if any
#' elements in `x` is negative, or if `x` is not a symmetric matrix or an
#' array of symmetric matrices.
check_dist <- function(x, name = "x") {
    if (any(is.na(x))) {
        stop("NA or NaN found in `", name, "`.", call. = FALSE)
    }

    if (any(x < 0)) {
        stop("invalid negative distance in `", name, "`.", call. = FALSE)
    }

    if (!is_numeric_scalar(x)) {
        if (is.array(x)) {
            if (!(length(dim(x)) %in% 2:3)) {
                stop("`", name, "` must be a matrix or 3-d array.",
                     call. = FALSE)
            }

            if (is.matrix(x)) {
                if (!isSymmetric.matrix(x))
                    stop("distance matrix `", name, "` is not symmetric.",
                         call. = FALSE)
            } else {
                for (i in 1:dim(x)[3])
                    if (!isSymmetric.matrix(x[, , i]))
                        stop("not all matrix slices in array `", name,
                             "` is symmetric.", call. = FALSE)
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
#' Check if `x` is a valid signed distance vector, matrix or array. It errors
#' if `x` in absolute value is not a symmetric matrix or an array of
#' symmetric matrices.
check_dist_sign <- function(x, name) {
    if (any(is.na(x))) {
        stop("NA or NaN found in `", name, "`.", call. = FALSE)
    }

        if (!is_numeric_scalar(x)) {
        if (is.array(x)) {
            if (!(length(dim(x)) %in% 2:3)) {
                stop("`", name, "` must be a matrix or 3-d array.",
                     call. = FALSE)
            }

            if (is.matrix(x)) {
                if (!isSymmetric.matrix(abs(x)))
                    stop("distance matrix `", name,
                         "` is not symmetric in absolute values.",
                         call. = FALSE)
            } else {
                for (i in 1:dim(x)[3])
                    if (!isSymmetric.matrix(abs(x)[, , i]))
                        stop("not all matrix slices in array `", name,
                             "`` is symmetric in absolute values.",
                             call. = FALSE)
            }
        }
    }
    return(invisible(NULL))
}

#' Check if valid input length
#'
#' @param x Scaler or vector
#' @param length Length of `x`.
#' @param name Name of `x`.
#'
#' @return `x`.
#' @keywords internal
#'
#' @details
#' Check if `x` has approprate length. If length of `x` is 1 then `x` is
#' replicated to match `length`. If length of `x` is neither 1 or `length` then
#' an error is signaled.
check_length <- function(x, length, name) {
    if (length(x) == 1) {
        x <- rep(x, length)
    } else {
        if (length(x) != length)
            stop("length of `", name, "` must be 1 or ", length, ".",
                 call. = FALSE)
    }
    return(x)
}

#' Check if valid input length
#'
#' @param x_ls List of scaler or vector
#' @param length_ls List of length of `x_ls`.
#' @param name Mame of `x_ls`.
#'
#' @return `x_ls`.
#' @keywords internal
#'
#' @details
#' Check if elements in `x_ls` have approprate length. If length of any elements
#' in `x_ls` is 1 then they are replicated to match `length`. If length of any
#' elements is neither 1 or `length` then an error is signaled.
check_length_ls <- function(x_ls, length, name) {

    for (k in 1:length(x_ls)) {
        if (length(x_ls[[k]]) == 1) {
            x_ls[[k]] <- rep(x_ls[[k]], length)
        } else {
            if (length(x_ls[[k]]) != length) {
                stop("length of `", name, "` must be 1 or ", length,
                     " for regime ", k, ".", call. = FALSE)
            }
        }
    }
    return(x_ls)
}
