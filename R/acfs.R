#' Generic function for calculating autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Mean of autocorrelations for each time lag.
#' @export
#' @family {functions related to the class}
acfs <- function(x, ...) {
    UseMethod("acfs")
}

#' Calculating mean autocorrelations for an `mcgf` object
#'
#' @name acfs.mcgf
#' @aliases `acfs<-`
#'
#' @param x An `mcgf` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param ... Additional parameters or attributes.
#'
#' @return Mean autocorrelations.attribute `acfs`.
#' @export
#'
#' @details
#' It computes mean autocorrelations for each time lag across locations. Use
#' [`add_acfs()`] to add `acfs` to `x`.
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' acfs(x = wind_mcgf, lag_max = 3)
#' @family {functions related to the class}
acfs.mcgf <- function(x, lag_max, ...) {
    acfs <- attr(x, "acfs")

    if (!is.null(acfs)) {
        return(acfs)
    } else {
        ccfs <- attr(x, "ccfs")
        if (!is.null(ccfs) && dim(ccfs)[3] != lag_max + 1)
            warning("`lag_max` must be the same as that in `ccfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.")

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.")

        data <- x
        n_var <- ncol(data)

        acf_data <- matrix(NA, nrow = lag_max + 1, ncol = n_var)

        for (i in 1:n_var) {
            acf_data[, i] <- stats::acf(data[, i], lag.max = lag_max,
                                        plot = FALSE, ...)$acf
        }

        acfs <- rowMeans(acf_data)
        names(acfs) <- paste0("lag", 0:lag_max)
        return(acfs)
    }
}

#' @rdname acfs.mcgf
#' @param value A Vector of mean of autocorrelations for time lags starting
#' from 0.
#' @export
`acfs<-` <- function(x, value) {

    attr(x, "acfs") <- value
    return(x)
}

#' Generic function for calculating and adding autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return `x` with mean of autocorrelations for each time lag.
#' @export
#' @family {functions related to the class}
add_acfs <- function(x, ...) {
    UseMethod("add_acfs")
}

#' Calculating and adding mean autocorrelations for an `mcgf` object
#'
#' @param x An `mcgf` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param ... Additional parameters or attributes.
#'
#' @return An `mcgf` object with newly added mean autocorrelations under
#' attribute `acfs`.
#' @export
#'
#' @details
#' It computes mean autocorrelations for each time lag across locations, and add
#' them to `x`.
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' wind_mcgf <- add_acfs(x = wind_mcgf, lag_max = 3)
#' print(wind_mcgf, "acfs")
#' @family {functions related to the class}
add_acfs.mcgf <- function(x, lag_max, ...) {

    acfs <- acfs.mcgf(x = x, lag_max = lag_max, ...)
    attr(x, "acfs") <- acfs
    return(x)
}
