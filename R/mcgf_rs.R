#' Create an mcgf_rs object
#'
#' @param x An mcgf object.
#' @param label A vector of regime labels. Its length must be the same as
#' the number rows in `data`.
#'
#' @keywords internal
#' @return An S3 object of class `mcgf_rs`. As it inherits and extends the
#' `mcgf` and then the`data.frame` class, all methods remain valid to the
#' `data` part of the object. Additional attributes may be assigned and
#' extracted.
new_mcgf_rs <- function(x, label) {

    data <- as.data.frame(x)
    label <- as.factor(label)

    return(structure(.Data = data, label = label,
                  class = c("mcgf_rs", "mcgf", "data.frame")))
}

#' Validate an mcgf_rs object
#'
#' @param An mcgf_rs object
#'
#' @keywords internal
#' @return An S3 object of class `mcgf_rs`.
#'
#' @details
#' It validates an `mcgf_rs` object by checking if `label` is of the matching
#' length.
validate_mcgf_rs <- function(x) {

    n <- nrow(x)
    label <- attr(x, "label", exact = TRUE)

    if (length(label) != n)
        stop("length of `label` must be the same as the number of rows in `x`.")

    return(x)
}

#' Create mcgf_rs object
#'
#' @inherit new_mcgf_rs return params return
#'
#' @param data Time series data set in space-wide format.
#' @param locations  A matrix of data.frame of 2D points, first column
#' longitude, second column latitude, both in decimal degrees. Required when
#' `dists` is not supplied.
#' @param dists List of signed distance matrices. Required when `locations` is
#' not supplied.
#' @param time Optional, a vector of equally spaced time stamps.
#'
#' @export
#'
#' @details
#' An `mcgf_rs` object extends the S3 classes `mcgf` and `data.frame`.
#'
#' For inputs, `data` must be in space-wide format where rows correspond to
#' different time stamps and columns refer to spatial locations. Supply either
#' `locations` or `dists`. `locations` is a matrix or data.frame of 2D points
#' with first column longitude and second column latitude. Both columns must be
#' in decimal degrees. Number of rows in `locations` must be the same as the
#' number of columns of `data`. `dists` must be a list of signed distance
#' matrices with names `h1`, `h2`, and `h`. If `h` is not given, it will be
#' calculated as the Euclidean distance of `h1` and `h2`. `time` is a vector of
#' equally spaced time stamps. If it is not supplied then `data` is assumed to
#' be temporally equally spaced. `label` must be a vector containing regime
#' labels, and its length must be the same as the number of rows in `x`.
#'
#' An `mcgf_rs` object extends the S3 classes `mcgf` and `data.frame`, all
#' methods remain valid to the `data` part of the object.
#'
#' @examples
#' data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
#' lon <- c(110, 120, 130)
#' lat <- c(50, 55, 60)
#' locations <- cbind(lon, lat)
#' label <- c(1, 1, 2, 2, 2)
#' obj <- mcgf_rs(data, locations = locations, label = label)
#' print(obj, "locations")
#' print(obj, "label")
#'
#' @family {functions related to the rs class}
mcgf_rs <- function(data, locations, dists, label, time) {

    if (length(unique(label)) == 1)
        message("only 1 unique class found in `label`, consider using `mcgf()`",
                "instead.")

    x_mcgf <- mcgf(data = data, locations = locations, dists = dists,
                   time = time)


    return(validate_mcgf_rs(new_mcgf_rs(x_mcgf, label)))
}

#' Check if an object is an `mcgf_rs` object..
#'
#' @name is.mcgf_rs
#' @param x An Object.
#'
#' @return `is.mcgf_rs` returns a logical valud; TRUE if `x` is of the `mcgf_rs`
#' class. `as.mcgf_rs` coerces an `mcgf` object to an `mcgf_rs` object by adding
#' regime labels.
#' @export
is.mcgf_rs <- function(x){
    inherits(x, "mcgf_rs")
}

#' @rdname is.mcgf_rs
#' @param label A vector of regime labels. Its length must be the same as
#' the number rows in `data`.
as.mcgf_rs <- function(x, label) {
    return(validate_mcgf_rs(new_mcgf_rs(x, label)))
}
