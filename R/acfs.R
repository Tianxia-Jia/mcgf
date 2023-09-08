#' Calculate regime-switching auto-correlation
#'
#' @param x A univariate numeric time series.
#' @param label A factor of regime labels.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param demean Logical. Should the covariances be about the sample means?
#'
#' @return Mean auto-correlations for each group in `label`.
#' @export
acf_rs <- function(x, label, lag_max, demean = TRUE) {

    stopifnot(length(x) == length(label))

    if (demean) {
        x <- x - mean(x)
    }

    n_reg <- length(unique(label))
    n_x <- length(x)
    lvs <- levels(label)

    lag_max <- ifelse(lag_max >= n_x, n_x - 1, lag_max)

    x <- stats::ts(x)
    acf_ls <- lapply(1:n_reg, function(x) {
        x <- numeric(lag_max + 1)
        x
    })

    for (u in 0:lag_max) {
        x_u <- stats::lag(x, -u)
        x_x_u <- (x * x_u)
        label_u <- label[(1 + u):n_x]

        for(k in 1:n_reg) {
            numer <- x_x_u[label_u == lvs[k]]
            if(length(numer) == 0) {
                acf_ls[[k]][u + 1] <- NA
            } else {
                acf_ls[[k]][u + 1] <- sum(numer)
            }
        }
    }

    acf_ls <- lapply(acf_ls, function(x) x / x[1])
    acf_ls <- lapply(acf_ls, function(x) {
        names(x) <- paste0("lag", 0:lag_max)
        x
    })
    names(acf_ls) <- paste0("Regime ", lvs)
    return(acf_ls)
}

#' Generic function for calculating autocorrelation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @details
#' Refer to [`acfs.mcgf()`] for more details.
#'
#' @export
acfs <- function(x, ...) {
    UseMethod("acfs")
}

#' Extract, calculate, or assign mean auto-correlations for an `mcgf` or
#' `mcgf_rs` object
#'
#' @name acfs.mcgf
#'
#' @param x An `mcgf` or `mcgf_rs` object.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param replace Logical; if TRUE, `acfs` are recalculated.
#' @param ... Additional parameters or attributes.
#'
#' @return [`acfs()`] returns (regime-switching) mean auto-correlations.
#' [`add_acfs()`] returns the same object with additional attributes of
#' (regime-switching) mean auto-correlations.
#' @export
#'
#' @details
#'
#' For `mcgf` objects, [`acfs()`] computes mean auto-correlations for each time
#' lag across locations. The output is a vector of acfs.
#'
#' For `mcgf` objects, [`acfs()`] computes regime-switching mean
#' auto-correlations for each time lag across locations. The output is a list of
#' vectors of acfs, where each vector in the list corresponds to the acfs for
#' a regime.
#'
#' [`acfs<-`] assigns `acfs` to `x`.
#'
#' [`add_acfs()`] adds `acfs` to `x`.
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' acfs(x = wind_mcgf, lag_max = 3)
#'
#' @family {functions related to the auto- and cross-correlations}
acfs.mcgf <- function(x, lag_max, replace = FALSE, ...) {
    acfs <- attr(x, "acfs", exact = TRUE)

    if (!is.null(acfs) && !replace) {
        return(acfs)
    } else {
        ccfs <- attr(x, "ccfs", exact = TRUE)
        if (!is.null(ccfs) && !is.mcgf_rs(x) && dim(ccfs)[3] != lag_max + 1)
            warning("`lag_max` must be the same as that in `ccfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.", call. = FALSE)

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.", call. = FALSE)

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
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf_rs(data = wind_sq, locations = wind_loc, time = time,
#' label = c(rep(1,3574), rep(2, 3000)))
#' acfs(x = wind_mcgf, lag_max = 3)
#' @export
acfs.mcgf_rs <- function(x, lag_max, replace = FALSE, ...) {
    acfs <- attr(x, "acfs", exact = TRUE)

    if (!is.null(acfs) && !replace) {
        return(acfs)
    } else {
        label <- attr(x, "label", exact = TRUE)
        ccfs <- attr(x, "ccfs", exact = TRUE)
        if (!is.null(ccfs) && dim(ccfs$ccfs)[3] != lag_max + 1)
            warning("`lag_max` must be the same as that in `ccfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.", call. = FALSE)

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.", call. = FALSE)

        data <- x
        n_var <- ncol(data)

        acf_data <- list()
        acf_data <- acf_rs(data[, 1], label = label, lag_max = lag_max)

        for (i in 2:n_var) {
            z <- acf_rs(data[, i], label = label, lag_max = lag_max)
            acf_data <- Map(cbind, acf_data, z)
        }

        acfs_rs <- lapply(acf_data, rowMeans)
        acfs <- acfs.mcgf(x, lag_max = lag_max)
        return(list(acfs = acfs, acfs_rs = acfs_rs))
    }
}

#' @rdname acfs.mcgf
#' @param value A Vector of mean of auto-correlations for time lags starting
#' from 0.
#' @export
`acfs<-` <- function(x, value) {

    attr(x, "acfs") <- value
    return(x)
}

#' @rdname acfs.mcgf
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' wind_mcgf <- add_acfs(x = wind_mcgf, lag_max = 3)
#' print(wind_mcgf, "acfs")
#'
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf_rs(data = wind_sq, locations = wind_loc, time = time,
#' label = c(rep(1,3574), rep(2, 3000)))
#' wind_mcgf <- add_acfs(x = wind_mcgf, lag_max = 3)
#' print(wind_mcgf, "acfs")
#' @export
add_acfs <- function(x, lag_max, ...) {

    acfs <- acfs(x = x, lag_max = lag_max, ...)
    attr(x, "acfs") <- acfs
    return(x)
}
