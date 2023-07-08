#' Generic function for calculating autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Mean of auto-correlations.
#' @export
#' @family {functions related to the class}
add_acfs <- function(x, ...) {
    UseMethod("add_acfs")
}

#' Calculating mean autocorrelations for an `mcgf` object
#'
#' @param x An `mcgf` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of mean autocorrelations.
#' @export
#'
#' @details
#' It computes mean autocorrelations for each time lag across locations, and add
#' them to `x`.
#'
#' @examples
#' data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
#' lon <- c(110, 120, 130)
#' lat <- c(50, 55, 60)
#' locations <- cbind(lon, lat)
#' obj <- mcgf(data = data, locations = locations)
#' obj <- add_acfs(x = obj, lag_max = 2)
#' print(obj, "acfs")
#' @family {functions related to the class}
add_acfs.mcgf <- function(x, lag_max, ...) {

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
    attr(x, "acfs") <- acfs
    return(x)
}
