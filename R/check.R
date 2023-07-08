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
#' @return NULL.
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
#' @return NULL.
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
#' @param name Name of `x_ls`.
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

#' Check if valid dists attribute for an `mcgf` object
#'
#' @param dists List of scaler or vector
#' @param n_var Scaler, number of variables.
#' @param name names of `dists`.
#'
#' @return `dists`.
#' @keywords internal
#'
#' @details
#' Check if `dists` is a valid dists attribute for an `mcgf` object. It errors
#' if 1) `dists` does not contain `h1` or `h2`, 2) if their dimensions do not
#' match, 3) if it contains elements other than `h`, `h1` or `h2`. `h` will be
#' computed if it is not given.
check_dists <- function(dists, n_var, names, name_dists = "dists") {

    if (!is.matrix(dists$h1))
        stop("`h1` in `", name_dists, "` must be a matrix.", call. = FALSE)
    if (!is.matrix(dists$h2))
        stop("`h2` in `", name_dists, "` must be a matrix.", call. = FALSE)

    if (any(dim(dists$h1) != c(n_var, n_var)))
        stop("`h1` in `", name_dists, "` must be a matrix of dimension ",
             n_var, " x ", n_var, ".", call. = FALSE)
    if (any(dim(dists$h2) != c(n_var, n_var)))
        stop("`h2` in `", name_dists, "` must be a matrix of dimension ",
             n_var, " x ", n_var, ".", call. = FALSE)

    check_dist_sign(dists$h1, name = "h1")
    check_dist_sign(dists$h2, name = "h2")

    if (is.null(rownames(dists$h1)))
        rownames(dists$h1) <- names

    if (is.null(colnames(dists$h1)))
        colnames(dists$h1) <- names

    if (is.null(rownames(dists$h2)))
        rownames(dists$h2) <- names

    if (is.null(colnames(dists$h2)))
        colnames(dists$h2) <- names

    if (!is.null(dists$h)) {
        if (!is.matrix(dists$h))
            stop("'h' in `", name_dists, "` must be a matrix.", call. = FALSE)

        if (any(dim(dists$h) != c(n_var, n_var)))
            stop("'h' in `", name_dists, "` must be a matrix of dimension ",
                 n_var, " x ", n_var, ".", call. = FALSE)

        check_dist(x = dists$h, "h")

        if (is.null(rownames(dists$h)))
            rownames(dists$h) <- names

        if (is.null(colnames(dists$h)))
            colnames(dists$h) <- names

    } else {
        dists$h <- sqrt(dists$h1 ^ 2 + dists$h2 ^ 2)
    }

    if (any(!names(dists) %in% c("h", "h1", "h2"))) {
        stop("invalid dists attribute in `", name_dists, "`.")
    }
    return(dists)
}
