#' Generic function for calculating distance matrices
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A list of signed distance matrices: `h` (Euclidean), `h1`
#' (horizontal), and `h2` (vertical) with the same dimensions.
#' @export
#' @family functions related to the class
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
#' @param return_grid Logical; used when `locations` in `x` are longitudes and
#' latitudes.
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
#'
#' # if locations are longitudes and latitudes
#' obj <- mcgf(data = data, locations = locations)
#' obj
#' dists(obj)
#' dists(obj) <- dists(obj)
#' obj
#'
#' # if locations are just coordinates in a 2D plane:
#' obj <- mcgf(data = data, locations = locations, longlat = FALSE)
#' obj
#'
#' # calculate distances
#' dists(obj)
#'
#' # add distances to the `mcgf` object
#' dists(obj) <- dists(obj)
#' obj
dists.mcgf <- function(x, return_grid = FALSE, ...) {
    if (return_grid && attr(x, "longlat", exact = TRUE)) {
        locations <- attr(x, "locations", exact = TRUE)
        origin <- attr(x, "origin", exact = TRUE)

        dists <- find_dists(locations,
            longlat = TRUE, origin = origin,
            return_grid = TRUE, ...
        )
    } else {
        dists <- attr(x, "dists", exact = TRUE)
    }

    if (is.null(dists)) {
        locations <- attr(x, "locations", exact = TRUE)
        longlat <- attr(x, "longlat", exact = TRUE)
        origin <- attr(x, "origin", exact = TRUE)

        if (return_grid) {
            dists <- find_dists(locations,
                longlat = longlat, origin = origin,
                return_grid = TRUE, ...
            )
        } else {
            dists <- find_dists(locations,
                longlat = longlat, origin = origin,
                ...
            )
        }
    }
    return(dists)
}

#' @rdname dists.mcgf
#' @param value List of signed distance matrices, outputted from [dists()].
#' @export
`dists<-` <- function(x, value) {
    dists <- check_dists(
        dists = value, n_var = ncol(x),
        names = colnames(x), name_dists = "value"
    )
    attr(x, "dists") <- value
    return(x)
}

#' Generate random distance matrices
#'
#' @param N Number of locations.
#' @param names Names of locations.
#' @param scale Scale of the distance matrices. Default is 100.
#'
#' @return List of signed distances.
#' @export
#'
#' @details
#' This function generates random distance matrices using `rnorm`. `scale`
#' controls the scale of the distance matrices.
#'
#' @examples
#' set.seed(123)
#' rdists(3)
#' rdists(3, scale = 1)
#' rdists(3, names = LETTERS[1:3])
rdists <- function(N, names, scale = 100) {
    if (!is_numeric_scalar(N)) {
        stop("`N` must be an integer.")
    }

    x <- stats::rnorm(N) * scale
    y <- stats::rnorm(N) * scale
    grid <- cbind(x, y)

    if (!missing(names)) {
        if (length(names) != N) {
            stop("Length of `names` must be `N`.")
        }
    } else {
        names <- paste0("loc", 1:N)
    }

    dists <- .find_dists(grid, names = names, longlat = F)
    return(dists)
}
