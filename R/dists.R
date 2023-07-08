#' Generic function for calculating distance matrices
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A list of signed distance matrices: `h` (Euclidean), `h1`
#' (horizontal), and `h2` (vertical).
#' @export
#' @family {functions related to the class}
dists <- function(x, ...) {
    UseMethod("dists")
}

#' Calculating distance matrices for an `mcgf` object
#'
#' @name dists.mcgf
#' @aliases `dists<-`
#'
#' @param x An `mcgf` object.
#' @param ... Additional parameters or attributes.
#'
#' @return A list of signed distance matrices: `h` (Euclidean), `h1`
#' (horizontal), and `h2` (vertical).
#' @export
#'
#' @details
#' If the `dists` attribute is available in `x`, it will be printed. Otherwise
#' `dists` will be calculated based on the `locations` attribute.
#'
#' @examples
#' data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
#' lon <- c(110, 120, 130)
#' lat <- c(50, 55, 60)
#' locations <- cbind(lon, lat)
#' obj <- mcgf(data = data, locations = locations)
#' obj
#' dists(obj)
#' dists(obj) <- dists(obj)
#' obj
#' @family {functions related to the class}
dists.mcgf <- function(x, ...) {

    dists <- attr(x, "dists")

    if (is.null(dists)) {
        locations <- attr(x, "locations")
        dists <- find_dists(locations, ...)
    }
    return(dists)
}

#' @rdname dists.mcgf
#' @param value List of signed distance matrices, outputted from [dists()].
#' @export
`dists<-` <- function(x, value) {

    dists <- check_dists(dists = value, n_var = ncol(x),
                         names = colnames(x), name_dists = "value")
    attr(x, "dists") <- value
    return(x)
}
