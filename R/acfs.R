#' Calculate regime-switching auto-correlation
#'
#' @param x A univariate numeric time series.
#' @param label A factor of regime labels.
#' @param lag_max Maximum lag at which to calculate the acf.
#' @param demean Logical. Should the covariances be about the sample means?
#'
#' @return Mean auto-correlations for each group in `label`.
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' label <- sample(1:2, 100, replace = TRUE)
#' acf_rs(x, label = factor(label), lag_max = 3)
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

        for (k in 1:n_reg) {
            numer <- x_x_u[label_u == lvs[k]]
            if (length(numer) == 0) {
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
#' @return A vector of mean auto-correlations for `mcgf` objects, or that plus
#' a list of regime-switching mean auto-correlations for `mcgf_rs` objects.
#'
#' @details
#' Refer to [`acfs.mcgf()`] and [`acfs.mcgf_rs()`] for more details.
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
#' For `mcgf_rs` objects, [`acfs()`] computes regime-switching mean
#' auto-correlations for each time lag across locations. The output is a list of
#' vectors of acfs, where each vector in the list corresponds to the acfs for
#' a regime.
#'
#' [`acfs<-`] assigns `acfs` to `x`.
#'
#' [`add_acfs()`] adds `acfs` to `x`.
#'
#' @examples
#' # Calculate acfs for 'sim1'
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#' acfs(sim1_mcgf, lag_max = 5)
#'
#' # Add acfs to 'sim1_mcgf'
#' sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
#' print(sim1_mcgf, "acfs")
#'
#' # Calculate acfs for 'sim2'
#' data(sim2)
#' sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#' acfs(sim2_mcgf, lag_max = 5)
#'
#' # Add acfs to 'sim2_mcgf'
#' sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
#' print(sim2_mcgf, "acfs")
#' @family functions related to calculating acfs and ccfs
acfs.mcgf <- function(x, lag_max, replace = FALSE, ...) {
    acfs <- attr(x, "acfs", exact = TRUE)

    if (!is.null(acfs) && !replace) {
        return(acfs)
    } else {
        ccfs <- attr(x, "ccfs", exact = TRUE)
        if (!is.null(ccfs) && !is.mcgf_rs(x) && dim(ccfs)[3] != lag_max + 1) {
            warning("`lag_max` must be the same as that in `ccfs`")
        }

        if (!is_numeric_scalar(lag_max)) {
            stop("`lag_max` must be numeric.", call. = FALSE)
        }

        if (lag_max < 0) {
            stop("`lag_max` must be a positive integer.", call. = FALSE)
        }

        data <- x
        n_var <- ncol(data)

        acf_data <- matrix(NA, nrow = lag_max + 1, ncol = n_var)

        for (i in 1:n_var) {
            acf_data[, i] <- stats::acf(data[, i],
                lag.max = lag_max,
                plot = FALSE, ...
            )$acf
        }

        acfs <- rowMeans(acf_data)
        names(acfs) <- paste0("lag", 0:lag_max)
        return(acfs)
    }
}

#' @rdname acfs.mcgf
#' @export
acfs.mcgf_rs <- function(x, lag_max, replace = FALSE, ...) {
    acfs <- attr(x, "acfs", exact = TRUE)

    if (!is.null(acfs) && !replace) {
        return(acfs)
    } else {
        label <- attr(x, "label", exact = TRUE)
        ccfs <- attr(x, "ccfs", exact = TRUE)
        if (!is.null(ccfs) && dim(ccfs$ccfs)[3] != lag_max + 1) {
            warning("`lag_max` must be the same as that in `ccfs`")
        }

        if (!is_numeric_scalar(lag_max)) {
            stop("`lag_max` must be numeric.", call. = FALSE)
        }

        if (lag_max < 0) {
            stop("`lag_max` must be a positive integer.", call. = FALSE)
        }

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
#' @export
add_acfs <- function(x, lag_max, ...) {
    acfs <- acfs(x = x, lag_max = lag_max, ...)
    attr(x, "acfs") <- acfs
    return(x)
}
