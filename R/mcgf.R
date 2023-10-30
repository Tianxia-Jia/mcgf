#' Create an mcgf object
#'
#' @param data Time series data set in space-wide format.
#' @param locations  A matrix of data.frame of 2D points, first column
#' longitude, second column latitude, both in decimal degrees. Required when
#' `dists` is not supplied.
#' @param dists List of signed distance matrices. Required when `locations` is
#' not supplied.
#' @param time Optional, a vector of equally spaced time stamps.
#'
#' @keywords internal
#' @return An S3 object of class `mcgf`. As it inherits and extends the
#' `data.frame` class, all methods remain valid to the `data` part of the
#' object. Additional attributes may be assigned and extracted.
new_mcgf <- function(data, locations, dists, time) {
    data <- as.data.frame(data)
    rownames(data) <- time

    if (!missing(dists)) {
        structure(.Data = data, dists = dists, class = c("mcgf", "data.frame"))
    } else {
        structure(
            .Data = data, locations = locations,
            class = c("mcgf", "data.frame")
        )
    }
}

#' Validate an mcgf object
#'
#' @param x An mcgf object.
#'
#' @keywords internal
#' @return An S3 object of class `mcgf`.
#'
#' @details
#' It validates an `mcgf` object by checking if `dists` contains valid
#' distance matrics/arrays.
validate_mcgf <- function(x) {
    data <- x
    n_var <- ncol(x)
    locations <- attr(x, "locations", exact = TRUE)
    dists <- attr(x, "dists", exact = TRUE)

    if (!is.null(dists)) {
        dists <- check_dists(
            dists = dists, n_var = n_var,
            names = colnames(data)
        )
        attr(x, "dists") <- dists
    }
    return(x)
}

#' Create mcgf object
#'
#' @inherit new_mcgf return params return
#'
#' @export
#'
#' @details
#' An `mcgf` object extends the S3 class `data.frame`.
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
#' be temporally equally spaced.
#'
#' An `mcgf` object extends the S3 class `data.frame`, all methods remain valid
#' to the `data` part of the object.
#'
#' @examples
#' data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
#' lon <- c(110, 120, 130)
#' lat <- c(50, 55, 60)
#' locations <- cbind(lon, lat)
#' obj <- mcgf(data, locations = locations)
#' print(obj, "locations")
mcgf <- function(data, locations, dists, time) {
    if (!is.data.frame(data) && !is.matrix(data)) {
        stop("`data` must be a matrix or data.frame.", call. = FALSE)
    }

    if (any(is.na(data))) {
        stop("`data` must not contain missing values.", call. = FALSE)
    }

    if (any(sapply(data, function(x) !is.numeric(x)))) {
        stop("non numeric values found in `data`.", call. = FALSE)
    }

    if (missing(locations) && missing(dists)) {
        stop("must provide either `locations` or `dists`.",
            call. = FALSE
        )
    }

    if (!missing(locations) && !missing(dists)) {
        stop("do not provide both `locations` and `dists`.", call. = FALSE)
    }

    name_var <- colnames(data)

    n_var <- NCOL(data)

    if (missing(time)) {
        message(
            "`time` is not provided, ",
            "assuming rows are equally spaced temporally."
        )
        time <- 1:NROW(data)
    }

    if (length(time) != NROW(data)) {
        stop("length of `time` must be the same as the number of rows of ",
            "`data`.",
            call. = FALSE
        )
    }

    diff_time <- diff(time)
    if (length(unique(diff_time)) != 1) {
        stop("`time` must be equally spaced.")
    }
    if (unique(diff_time) < 0) {
        stop("`time` must be in ascending order.")
    }

    if (!missing(locations)) {
        if (ncol(data) != nrow(locations)) {
            stop("number of columns of `data` must be the same as the ",
                "number of rows of `locations`",
                call. = FALSE
            )
        }

        if (!is.null(rownames(locations))) {
            locations <- locations[match(name_var, rownames(locations)), ]
        }

        if (any(colnames(data) != rownames(locations))) {
            stop("row names of `locations` are not the same as the column ",
                "names of `data`.",
                call. = FALSE
            )
            rownames(locations) <- colnames(data)
        }
        return(validate_mcgf(new_mcgf(
            data = data, locations = locations,
            time = time
        )))
    } else {
        if (!is.list(dists)) {
            stop("`dists` must be a list.", call. = FALSE)
        }
        if (any(!c("h1", "h2") %in% names(dists))) {
            stop("`dists` must contain 'h1' and 'h2',", call. = FALSE)
        }

        return(validate_mcgf(new_mcgf(data = data, dists = dists, time = time)))
    }
}

#' Check if an object is an `mcgf` object.
#'
#' @param x An Object.
#'
#' @return Logical; TRUE if `x` is of the `mcgf` class
#' @export
#' @examples
#' data(sim1)
#' is.mcgf(sim1)
#'
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#' is.mcgf(sim1_mcgf)
is.mcgf <- function(x) {
    inherits(x, "mcgf")
}
