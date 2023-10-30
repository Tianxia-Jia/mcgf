#' Generic function for calculating cross-correlation
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return An array of cross-correlations for `mcgf` objects, or that plus
#' a list of regime-switching cross-correlations for `mcgf_rs` objects.
#'
#' @details
#' Refer to [`ccfs.mcgf()`] and [`ccfs.mcgf_rs()`] for more details.
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
#' @param replace Logical; if TRUE, `acfs` are recalculated.
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
#' matrices where each array in the list corresponds to the cross-correlation
#' for a regime.
#'
#' [`ccfs<-`] assigns `ccfs` to `x`.
#'
#' [`add_ccfs()`] adds `ccfs` and `sds` to `x`.
#'
#' @export
#' @examples
#' # Calculate ccfs for 'sim1'
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#' ccfs(sim1_mcgf, lag_max = 5)
#'
#' # To use multiple cores, use the `ncores` argument
#' ccfs(sim1_mcgf, lag_max = 5, ncores = 2)
#'
#' # Add ccfs and sds to 'sim1_mcgf'
#' sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = 5)
#' print(sim1_mcgf, "ccfs")
#' print(sim1_mcgf, "sds")
#'
#' # Calculate ccfs for 'sim2'
#' data(sim2)
#' sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#' ccfs(sim2_mcgf, lag_max = 5)
#'
#' # Add ccfs and sds to 'sim2_mcgf'
#' sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)
#' print(sim2_mcgf, "ccfs")
#' print(sim2_mcgf, "sds")
#' @family functions related to acfs and ccfs
ccfs.mcgf <- function(x, lag_max, ncores = 1, replace = FALSE, ...) {
    ccfs <- attr(x, "ccfs", exact = TRUE)

    if (!is.null(ccfs) && !replace) {
        return(ccfs)
    } else {
        acfs <- attr(x, "acfs", exact = TRUE)
        if (!is.null(acfs) && !is.mcgf_rs(x) && length(acfs) != lag_max + 1) {
            warning("`lag_max` must be the same as that in `acfs`")
        }

        if (!is_numeric_scalar(lag_max)) {
            stop("`lag_max` must be numeric.")
        }

        if (lag_max < 0) {
            stop("`lag_max` must be a positive integer.")
        }

        data <- x
        n_var <- ncol(data)

        if ((n_var > 30 && lag_max > 10) || n_var > 50) {
            if (ncores == 1) {
                message(
                    "Large dataset, this may take a while. Set `ncores` > 1 to",
                    " speed up."
                )
            }
        }

        ccfs <- array(
            NA,
            dim = c(n_var, n_var, lag_max + 1),
            dimnames = list(
                colnames(data),
                colnames(data),
                paste0("lag", 0:lag_max)
            )
        )

        if (ncores == 1) {
            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs[i, j, ] <- stats::ccf(data[, i],
                        data[, j],
                        lag.max = lag_max,
                        plot = F
                    )$acf[-c(1:lag_max)]
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
            if (ncores > parallel::detectCores()) {
                ncores <- parallel::detectCores()
            }

            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)

            `%dopar%` <- foreach::`%dopar%`

            ccfs_ls <- foreach::foreach(i = 1:n_var) %dopar% {
                ccfs_i <- vector("list", n_var)
                for (j in 1:n_var) {
                    ccfs_i[[j]] <- stats::ccf(data[, i],
                        data[, j],
                        lag.max = lag_max,
                        plot = F
                    )$acf[-c(1:lag_max)]
                }
                ccfs_i
            }
            parallel::stopCluster(cl)

            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs[i, j, ] <- ccfs_ls[[i]][[j]]
                }
            }
        }
        return(ccfs)
    }
}

#' @rdname ccfs.mcgf
#' @export
ccfs.mcgf_rs <- function(x, lag_max, ncores = 1, replace = FALSE, ...) {
    ccfs <- attr(x, "ccfs", exact = TRUE)

    if (!is.null(ccfs) && !replace) {
        return(ccfs)
    } else {
        label <- attr(x, "label", exact = TRUE)
        acfs <- attr(x, "acfs", exact = TRUE)
        if (!is.null(acfs) && length(acfs$acfs) != lag_max + 1) {
            warning("`lag_max` must be the same as that in `acfs`")
        }

        if (!is_numeric_scalar(lag_max)) {
            stop("`lag_max` must be numeric.")
        }

        if (lag_max < 0) {
            stop("`lag_max` must be a positive integer.")
        }

        data <- x
        n_var <- ncol(data)
        n_regime <- length(levels(label))

        if ((n_var > 30 && lag_max > 10) || n_var > 50) {
            if (ncores == 1) {
                message(
                    "Large dataset, this may take a while. Set `ncores` > 1 to",
                    "speed up.\n"
                )
            }
        }

        ccfs_rs <- array(
            NA,
            dim = c(n_var, n_var, lag_max + 1),
            dimnames = list(
                colnames(data),
                colnames(data),
                paste0("lag", 0:lag_max)
            )
        )
        ccfs_rs <- rep(list(ccfs_rs), n_regime)
        names(ccfs_rs) <- paste0("Regime ", levels(label))

        if (ncores == 1) {
            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    ccfs_i_j <- ccf_rs(data[, i],
                        data[, j],
                        label = label,
                        lag_max = lag_max
                    )
                    for (k in 1:n_regime) {
                        ccfs_rs[[k]][i, j, ] <- ccfs_i_j[[k]][-c(1:lag_max)]
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
            if (ncores > parallel::detectCores()) {
                ncores <- parallel::detectCores()
            }

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
                        lag_max = lag_max
                    )
                }
                ccfs_i
            }
            parallel::stopCluster(cl)

            for (i in 1:n_var) {
                for (j in 1:n_var) {
                    for (k in 1:n_regime) {
                        ccfs_rs[[k]][i, j, ] <-
                            ccfs_ls[[i]][[j]][[k]][-c(1:lag_max)]
                    }
                }
            }
        }

        ccfs <- ccfs.mcgf(x, lag_max = lag_max, ncores = ncores)
        return(list(ccfs = ccfs, ccfs_rs = ccfs_rs))
    }
}

#' Calculate regime-switching cross-correlation
#'
#' @param x,y A univariate numeric time series.
#' @param label A factor of regime labels.
#' @param lag_max Maximum lag at which to calculate the ccf.
#'
#' @return Cross-correlations for each group in `label`.
#' @export
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100)
#' label <- sample(1:2, 100, replace = TRUE)
#' ccf_rs(x, y, label = factor(label), lag_max = 3)
ccf_rs <- function(x, y, label, lag_max) {
    stopifnot(length(x) == length(y))
    stopifnot(length(y) == length(label))

    x <- x - mean(x)
    y <- y - mean(y)

    n_reg <- length(unique(label))
    n_x <- length(x)
    lvs <- levels(label)

    lag_max <- ifelse(lag_max >= n_x, n_x - 1, lag_max)

    x <- stats::ts(x)
    y <- stats::ts(y)
    ccf_ls <- lapply(1:n_reg, function(x) {
        x <- numeric(2 * lag_max + 1)
        x
    })

    for (k in 1:n_reg) {
        x_k <- x[label == lvs[k]]
        y_k <- y[label == lvs[k]]
        denom <- sqrt(sum(x_k^2) * sum(y_k^2))

        for (u in 0:lag_max) {
            y.u <- stats::lag(y, -u)
            x_x_u <- (x * y.u)
            label_u <- label[(1 + u):n_x]
            numer <- x_x_u[label_u == lvs[k]]

            if (length(numer) == 0) {
                ccf_ls[[k]][lag_max + u + 1] <- NA
            } else {
                ccf_ls[[k]][lag_max + u + 1] <- sum(numer) / denom
            }
        }

        for (u in 1:lag_max) {
            y.u <- stats::lag(y, u)

            x_x_u <- (x * y.u)
            label_u <- label[(1 + u):n_x]
            numer <- x_x_u[label_u == lvs[k]]

            if (length(numer) == 0) {
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
    names(ccf_ls) <- paste0("Regime ", lvs)
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
