#' Generic function for calculating cross-correlation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @details
#' Refer to [`ccfs.mcgf()`] for more details.
#'
#' @export
ccfs <- function(x, ...) {
    UseMethod("ccfs")
}

#' Extract, calculate, or assign cross-correlations for an `mcgf` or `mcgf_rs`
#' object
#'
#' @name ccfs.mcgf
#'
#' @param x An `mcgf` or `mcgf_rs` object.
#' @param lag_max Maximum lag at which to calculate the ccfs.
#' @param ncores Number of cpu cores used for computing. The `doParallel`
#' package is required when `ncores` > 1.
#' @param ... Additional parameters or attributes. Not in use.
#'
#' @return [`ccfs()`] returns (regime-switching) cross-correlations.
#' [`add_ccfs()`] returns the same object with additional attributes of
#' (regime-switching) cross-correlations and (regime-switching) empirical
#' standard deviations.
#'
#' @details
#' For `mcgf` objects, [`ccfs()`] computes cross-correlations for each time
#' lag. The output is an array of matrices where each matrix corresponds to the
#' cross-correlation for a time lag.
#'
#' For `mcgf_rs` objects, [`ccfs()`] computes regime-switching
#' cross-correlations for each time lag. The output is a list of array of
#' matrices where each array in the list corresponds to the ccfs for a regime.
#'
#' [`ccfs<-`] assigns `ccfs` to `x`.
#'
#' [`add_ccfs()`] adds `ccfs` and `sds` to `x`.
#'
#' @export
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' ccfs(x = wind_mcgf, lag_max = 3)
#' # ccfs(x = wind_mcgf, lag_max = 3, ncores = 10)
#'
#' @family {functions related to the auto- and cross-correlations}
ccfs.mcgf <- function(x, lag_max, ncores = 1,  ...) {

    ccfs <- attr(x, "ccfs", exact = TRUE)

    if (!is.null(ccfs)) {
        return(ccfs)
    } else {
        acfs <- attr(x, "acfs", exact = TRUE)
        if (!is.null(acfs) && length(acfs) != lag_max + 1)
            warning("`lag_max` must be the same as that in `acfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.")

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.")

        data <- x
        n_var <- ncol(data)

        if ((n_var > 30 && lag_max > 10) || n_var > 50)
            cat("Large dataset, this may take a while. Set `ncores` > 1 to",
                "speed up.\n")

        ccfs <- array(
            NA,
            dim = c(n_var, n_var, lag_max + 1),
            dimnames = list(colnames(data),
                            colnames(data),
                            paste0("lag", 0:lag_max))
        )

        if (ncores == 1) {
            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs[i, j,] <- stats::ccf(data[, i],
                                              data[, j],
                                              lag.max = lag_max,
                                              plot = F)$acf[-c(1:lag_max)]
                }
            }
        } else {
            if (!requireNamespace("doParallel", quietly = TRUE)) {
                stop("The `doParallel` package is required when `ncores` > 1.")
            }
            if (!requireNamespace("parallel", quietly = TRUE)) {
                stop("The `parallel` package is required when `ncores` > 1.")
            }
            if (!requireNamespace("foreach", quietly = TRUE)) {
                stop("The `foreach` package is required when `ncores` > 1.")
            }
            if (ncores > parallel::detectCores())
                ncores <- parallel::detectCores()

            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)

            `%dopar%` <- foreach::`%dopar%`

            ccfs_ls <- foreach::foreach(i = 1:n_var) %dopar% {
                ccfs_i <- vector("list", n_var)
                for (j in 1:n_var) {
                    ccfs_i[[j]] <- stats::ccf(data[, i],
                                              data[, j],
                                              lag.max = lag_max,
                                              plot = F)$acf[-c(1:lag_max)]
                }
                ccfs_i
            }
            parallel::stopCluster(cl)

            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs[i, j,] <- ccfs_ls[[i]][[j]]
                }
            }
        }
        return(ccfs)
    }
}

#' @rdname ccfs.mcgf
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf_rs(data = wind_sq, locations = wind_loc, time = time,
#' label = c(rep(1,3574), rep(2, 3000)))
#' ccfs(x = wind_mcgf, lag_max = 3)
#' # ccfs(x = wind_mcgf, lag_max = 3, ncores = 10)
#' @export
ccfs.mcgf_rs <- function(x, lag_max, ncores = 1, ...) {
    ccfs <- attr(x, "ccfs", exact = TRUE)

    if (!is.null(ccfs)) {
        return(ccfs)
    } else {
        label <- attr(x, "label", exact = TRUE)
        acfs <- attr(x, "acfs", exact = TRUE)
        if (!is.null(acfs) && length(acfs) != lag_max + 1)
            warning("`lag_max` must be the same as that in `acfs`")

        if (!is_numeric_scalar(lag_max))
            stop("`lag_max` must be numeric.")

        if (lag_max < 0)
            stop("`lag_max` must be a positive integer.")

        data <- x
        n_var <- ncol(data)
        n_regime <- length(levels(label))

        if ((n_var > 30 && lag_max > 10) || n_var > 50)
            cat("Large dataset, this may take a while. Set `ncores` > 1 to",
                "speed up.\n")

        ccfs <- array(
            NA,
            dim = c(n_var, n_var, lag_max + 1),
            dimnames = list(colnames(data),
                            colnames(data),
                            paste0("lag", 0:lag_max))
        )
        ccfs <- rep(list(ccfs), n_regime)
        names(ccfs) <- paste0("Regime ", levels(label))

        if (ncores == 1) {
            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs_i_j <- ccf_rs(data[, i],
                                       data[, j],
                                       label = label,
                                       lag_max = lag_max)
                    for (k in 1:n_regime) {
                        ccfs[[k]][i, j, ] <- ccfs_i_j[[k]][-c(1:lag_max)]
                    }
                }
            }
        } else {
            if (!requireNamespace("doParallel", quietly = TRUE)) {
                stop("The `doParallel` package is required when `ncores` > 1.")
            }
            if (!requireNamespace("parallel", quietly = TRUE)) {
                stop("The `parallel` package is required when `ncores` > 1.")
            }
            if (!requireNamespace("foreach", quietly = TRUE)) {
                stop("The `foreach` package is required when `ncores` > 1.")
            }
            if (ncores > parallel::detectCores())
                ncores <- parallel::detectCores()

            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)

            `%dopar%` <- foreach::`%dopar%`
            parallel::clusterExport(cl, "ccf_rs")

            ccfs_ls <- foreach::foreach(i = 1:n_var) %dopar% {
                ccfs_i <- vector("list", n_var)
                for (j in 1:n_var) {
                    ccfs_i[[j]] <- ccf_rs(data[, i],
                                          data[, j],
                                          label = label,
                                          lag_max = lag_max)
                }
                ccfs_i
            }
            parallel::stopCluster(cl)

            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    for (k in 1:n_regime) {
                        ccfs[[k]][i, j,] <-
                            ccfs_ls[[i]][[j]][[k]][-c(1:lag_max)]
                    }
                }
            }
        }
        return(ccfs)
    }
}

#' Calculate regime-switching cross-correlation
#'
#' @param x,y A univariate numeric time series.
#' @param label A vector of regime labels.
#' @param lag_max Maximum lag at which to calculate the ccf.
#'
#' @return Cross-correlations for each group in `label`.
#' @export
ccf_rs <- function(x, y, label, lag_max) {

    stopifnot(length(x) == length(y))
    stopifnot(length(y) == length(label))

    x <- x - mean(x)
    y <- y - mean(y)

    n_reg <- length(unique(label))
    n_x <- length(x)

    lag_max <- ifelse(lag_max >= n_x, n_x - 1, lag_max)

    x <- stats::ts(x)
    y <- stats::ts(y)
    ccf_ls <- lapply(1:n_reg, function(x) {
        x <- numeric(2 * lag_max + 1)
        x
    })

    for(k in 1:n_reg) {

        x_k <- x[label == k]
        y_k <- y[label == k]
        denom <- sqrt(sum(x_k ^ 2) * sum(y_k ^ 2))

        for (u in 0:lag_max) {
            y.u <- stats::lag(y, -u)
            x_x_u <- (x * y.u)
            label_u <- label[(1 + u):n_x]
            numer <- x_x_u[label_u == k]

            if(length(numer) == 0) {
                ccf_ls[[k]][lag_max + u + 1] <- NA
            } else {
                ccf_ls[[k]][lag_max + u + 1] <- sum(numer) / denom
            }
        }

        for (u in 1:lag_max) {
            y.u <- stats::lag(y, u)

            x_x_u <- (x * y.u)
            label_u <- label[(1 + u):n_x]
            numer <- x_x_u[label_u == k]

            if(length(numer) == 0) {
                ccf_ls[[k]][lag_max + 1 - u] <- NA
            } else {
                ccf_ls[[k]][lag_max + 1 - u] <- sum(numer) / denom
            }
        }
    }

    ccf_ls <- lapply(ccf_ls, function(x) {
        names(x) <- paste0(-lag_max:lag_max)
        x
    })
    return(ccf_ls)
}


#' @rdname ccfs.mcgf
#' @param value Cross-correlations.
#' @export
`ccfs<-` <- function(x, value) {

    attr(x, "ccfs") <- value
    return(x)
}

#' @rdname ccfs.mcgf
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' wind_mcgf <- add_ccfs(x = wind_mcgf, lag_max = 3)
#' print(wind_mcgf, "ccfs")
#' print(wind_mcgf, "sds")
#' @export
add_ccfs <- function(x, lag_max, ncores = 1, ...) {

    ccfs <- ccfs(x = x, lag_max = lag_max, ncores = ncores, ...)
    attr(x, "ccfs") <- ccfs

    sds <- attr(x, "sds", exact = TRUE)
    if (is.null(sds)) {
        attr(x, "sds") <- sds(x)
    }
    return(x)
}
