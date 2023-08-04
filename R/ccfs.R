#' Generic function for calculating cross-autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Cross-correlations.
#' @export
#'
#' @family {functions related to the auto- and cross-correlations}
ccfs <- function(x, ...) {
    UseMethod("ccfs")
}

#' Calculating cross-correlations for an `mcgf` object
#'
#' @name ccfs.mcgf
#' @aliases `ccfs<-`
#'
#' @param x An `mcgf` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param ... Additional parameters or attributes.
#'
#' @return  Cross-correlations.
#' @export
#'
#' @details
#' It computes cross-correlations for each time lag. The output is an array of
#' matrices where each matrix corresponds to the cross-correlation for a time
#' lag. Use [`add_ccfs()`] to add `ccfs` and `sds` to `x`.
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' ccfs(x = wind_mcgf, lag_max = 3)
#'
#' @family {functions related to the auto- and cross-correlations}
ccfs.mcgf <- function(x, lag_max, ...) {

    ccfs <- attr(x, "ccfs")

    if (!is.null(ccfs)) {
        return(ccfs)
    } else {
        acfs <- attr(x, "acfs")
        if (!is.null(acfs) && length(acfs) != lag_max + 1)
            warning("`lag_max` must be the same as that in `acfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.")

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.")

        data <- x
        n_var <- ncol(data)

        if (nrow(data) > 1000 && n_var > 30)
            cat("Large dataset, this may take a while.\n")

        ccfs <- array(
            NA,
            dim = c(n_var, n_var, lag_max + 1),
            dimnames = list(colnames(data),
                            colnames(data),
                            paste0("lag", 0:lag_max))
        )

        for (i in 1:n_var) {
            for (j in 1:n_var) {
                ccfs[i, j, ] <- stats::ccf(data[, i],
                                           data[, j],
                                           lag.max = lag_max,
                                           plot = F)$acf[-c(1:lag_max)]
            }
        }
        return(ccfs)
    }
}

#' @rdname ccfs.mcgf
#' @param value Cross-correlations.
#' @export
`ccfs<-` <- function(x, value) {

    attr(x, "ccfs") <- value
    return(x)
}

#' Generic function for calculating cross-autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return `x` with cross-correlations and standard deviations.
#' @export
#'
#' @family {functions related to the auto- and cross-correlations}
add_ccfs <- function(x, ...) {
    UseMethod("add_ccfs")
}

#' Calculating cross-correlations for an `mcgf` object
#'
#' @name add_ccfs.mcgf
#' @aliases `sds<-`
#'
#' @param x An `mcgf` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param ... Additional parameters or attributes.
#'
#' @return  An `mcgf` object with newly added cross-correlations and a vector
#' of standard deviations under attributes `ccfs` and `sds` respectively.
#' @export
#'
#' @details
#' It computes cross-correlations for each time lag, and add them to `x`.
#' The output is an array of matrices where each matrix corresponds to the
#' cross-correlation for a time lag, and a vector of empirical standard
#' deviations. If other standard deviation functions are used, they can be
#' replaced via [sds<-].
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' wind_mcgf <- add_ccfs(x = wind_mcgf, lag_max = 3)
#' print(wind_mcgf, "ccfs")
#' print(wind_mcgf, "sds")
#'
#' @family {functions related to the auto- and cross-correlations}
add_ccfs.mcgf <- function(x, lag_max, ...) {

    if (!is_numeric_scalar(lag_max))
        stop("`lag_max` must be numeric.")

    if (lag_max < 0)
        stop("`lag_max` must be a positive integer.")

    data <- x
    n_var <- ncol(data)

    if (nrow(data) > 1000 && n_var > 30)
        cat("Large dataset, this may take a while.\n")

    ccfs <- array(
        NA,
        dim = c(n_var, n_var, lag_max + 1),
        dimnames = list(colnames(data),
                        colnames(data),
                        paste0("lag", 0:lag_max))
    )

    for (i in 1:n_var) {
        for (j in 1:n_var) {
            ccfs[i, j,] <- stats::ccf(data[, i],
                                      data[, j],
                                      lag.max = lag_max,
                                      plot = F)$acf[-c(1:lag_max)]
        }
    }

    sds <- apply(data, 2, stats::sd)

    attr(x, "ccfs") <- ccfs
    attr(x, "sds") <- sds
    return(x)
}
