#' Generic function for computing kriging forecasts for new locations
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Kriging results of `x`.
#' @export
#'
#' @details
#' Refer to [`krige_new.mcgf()`] and [`krige_new.mcgf_rs()`] for more details.
krige_new <- function(x, ...) {
    UseMethod("krige_new")
}

#' Obtain kriging forecasts for new locations for an `mcgf` object.
#'
#' @param x An `mcgf` object.
#' @param newdata A data.frame with the same column names as `x`. If `newdata`
#' is missing the forecasts at the original data points are returned.
#' @param locations_new A matrix of data.frame of 2D points of new locations,
#' first column longitude, second column latitude, both in decimal degrees.
#' Supply only if `x` contains `locations`. Required when `dists_new` is not
#' supplied.
#' @param dists_new List of signed distance matrices (vectors) with names `h`,
#' `h1`, and 'h2' for all locations, with new locations in the end. Each
#' matrix must have the same number of columns. Required when `locations_new` is
#' not supplied.
#' @param newdata_new Optional; a data.frame with the same number of rows as
#' `newdata`. It contains the data of the new locations.
#' @param sds_new The standard deviations of the new locations. Default is 1.
#' @param model Which model to use. One of `all` or `base`.
#' @param interval Logical; if TRUE, prediction intervals are computed.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for
#' the intervals (if any) to be calculated. Used when `interval = TRUE`
#' @param ... Additional arguments. Not in use.
#'
#' @return A list of kriging forecasts (and prediction intervals) for all
#' locations.
#' @export
#'
#' @details
#' It produces simple kriging forecasts for a zero-mean mcgf for new locations
#' given their coordinates or relative distances. It supports kriging for the
#' `base` model and the `all` model which is the general stationary model with
#' the base and Lagrangian model from `x`.
#'
#' Users can either supply the coordinates via `locations_new`, or a list of
#' distance for all locations via `dists_new`, with new locations at the
#' end. `dists_new` will be used to calculate the new covariance matrices.
#' When `locations_new` is used, make sure `x` contains the attribute
#' `locations` of the coordinates of the old locations. When `dists_new` is
#' used, it should be a list of signed distance matrices of the same dimension,
#' where each row corresponds to the relative distances between a new location
#' and old locations in the same order as they appear in `x`.
#'
#' If data for the new locations are available, use `newdata_new` to include
#' them and they will be used to calculate the kriging forecasts for the new
#' locations; otherwise only data of the old locations will be used via
#' `newdata`.
#'
#' When `interval = TRUE`, confidence interval for each forecasts and each
#' horizon is given. Note that it does not compute confidence regions.
#'
#' @examples
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, locations = sim1$locations)
#' sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
#' sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = 5)
#'
#' # Fit a separable model and store it to 'sim1_mcgf'
#' fit_sep <- fit_base(
#'     sim1_mcgf,
#'     model = "sep",
#'     lag = 5,
#'     par_init = c(
#'         c = 0.001,
#'         gamma = 0.5,
#'         a = 0.3,
#'         alpha = 0.5
#'     ),
#'     par_fixed = c(nugget = 0)
#' )
#' sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep)
#'
#' # Fit a Lagrangian model
#' fit_lagr <- fit_lagr(
#'     sim1_mcgf,
#'     model = "lagr_tri",
#'     par_init = c(v1 = 300, v2 = 300, lambda = 0.15),
#'     par_fixed = c(k = 2)
#' )
#'
#' # Store the fitted Lagrangian model to 'sim1_mcgf'
#' sim1_mcgf <- add_lagr(sim1_mcgf, fit_lagr = fit_lagr)
#'
#' # Calculate the simple kriging predictions and intervals for all locations
#' locations_new <- rbind(c(-110, 55), c(-109, 54))
#' sim1_krige <- krige_new(sim1_mcgf,
#'     locations_new = locations_new,
#'     interval = TRUE
#' )
#' @family functions on fitting an mcgf
krige_new.mcgf <- function(x, newdata = NULL, locations_new = NULL,
                           dists_new = NULL, newdata_new = NULL,
                           sds_new = 1, model = c("all", "base"),
                           interval = FALSE, level = 0.95, ...) {
    model <- match.arg(model)
    dots <- list(...)

    no_newdata <- ifelse(is.null(newdata), TRUE, FALSE)
    no_newdata_new <- ifelse(is.null(newdata_new), TRUE, FALSE)

    no_locations_new <- ifelse(is.null(locations_new), TRUE, FALSE)
    no_dists_new <- ifelse(is.null(dists_new), TRUE, FALSE)

    if (model == "base") {
        if (is.null(attr(x, "base", exact = T))) {
            stop("Base model missing from `x`.", call. = FALSE)
        }
    } else {
        if (is.null(attr(x, "lagr", exact = T))) {
            stop("Lagrangian model missing from `x`.", call. = FALSE)
        }
    }

    if (no_locations_new && no_dists_new) {
        stop("must provide either `locations_new` or `dists_new`.",
            call. = FALSE
        )
    }

    if (!no_locations_new && !no_dists_new) {
        stop("do not provide both `locations_new` or `dists_new`.",
            call. = FALSE
        )
    }

    if (!no_locations_new) {
        locations <- attr(x, "locations", exact = TRUE)

        if (is.vector(locations_new)) {
            if (length(locations_new) != 2) {
                stop("incorrect dimension for `locations_new`", call. = FALSE)
            }
            locations_new <- matrix(locations_new, nrow = 1)
        }

        if (min(locations_new[, 1]) < min(locations[, 1]) ||
            max(locations_new[, 1]) > max(locations[, 1])) {
            message("Longitude for new locations is outside of the grid.")
        }
        if (min(locations_new[, 2]) < min(locations[, 2]) ||
            max(locations_new[, 2]) > max(locations[, 2])) {
            message("Latitude for new locations is outside of the grid.")
        }

        if (is.null(locations)) {
            stop("`locations` not found in 'x', supply 'dists_new' instead.",
                call. = FALSE
            )
        }

        n_var_new <- nrow(locations_new)
        names_new <- rownames(locations_new)

        if (!no_newdata_new) {
            if (ncol(newdata_new) != n_var_new) {
                stop("number of columns of `newdata_new` must be the same as ",
                    "the number of rows of `locations`",
                    call. = FALSE
                )
            }
        }

        longlat <- attr(x, "longlat", exact = TRUE)
        origin <- attr(x, "origin", exact = TRUE)
        dists_new <- find_dists_new(locations, locations_new[, 1:2],
            longlat = longlat, origin = origin
        )
    }

    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_max <- lag + horizon - 1
    n_block <- lag_max + 1
    n_var <- ncol(dists(x)$h)

    if (!no_dists_new) {
        if (!is.list(dists_new)) {
            stop("`dists_new` must be a list.", call. = FALSE)
        }
        if (any(!c("h", "h1", "h2") %in% names(dists_new))) {
            stop("`dists_new` must contain 'h', 'h1' and 'h2',", call. = FALSE)
        }
        if (any(dim(dists_new$h) != dim(dists_new$h1))) {
            stop("dimensions differ for 'h' and 'h1' in `dists_new`.",
                call. = FALSE
            )
        }
        if (any(dim(dists_new$h1) != dim(dists_new$h2))) {
            stop("dimensions differ for 'h1' and 'h2' in `dists_new`.",
                call. = FALSE
            )
        }
        if (is.array(dists_new$h1) && !is.matrix(dists_new$h1)) {
            if (length(dim(dists_new$h1)) != 3) {
                stop("'h1' and 'h2' must be matrices or 3D arrays in ",
                    "`dists_new`.",
                    call. = FALSE
                )
            }
        }
        n_var_new <- nrow(dists_new$h) - n_var
        names_new <- rownames(dists_new$h)[(n_var + 1):nrow(dists_new$h)]
    }

    if (is.null(names_new)) {
        names_new <- paste0("New_", seq_len(n_var_new))
    }

    if (length(sds_new) == 1) {
        sds_new <- rep(sds_new, n_var_new)
    }

    if (!no_newdata) {
        if (NCOL(newdata) != ncol(x)) {
            stop("unmatching number of columns for `newdata`.", call. = FALSE)
        }
        if (NROW(newdata) < lag) {
            stop("number of rows in `newdata` must be higher than `lag` ",
                lag, ".",
                call. = FALSE
            )
        }
    } else {
        newdata <- as.matrix(x)
    }

    if (!no_newdata_new) {
        if (ncol(newdata_new) != n_var_new) {
            stop("number of columns of `newdata_new` does not math with ",
                "`dists_new` or `locations_new`",
                call. = FALSE
            )
        }
        if (no_newdata) {
            if (nrow(newdata_new) != nrow(x)) {
                stop("number of rows do not match for `newdata_new` and `x`.",
                    call. = FALSE
                )
            }
        } else {
            if (nrow(newdata_new) != nrow(newdata)) {
                stop("number of rows do not match for `newdata_new` and ",
                    "`newdata`.",
                    call. = FALSE
                )
            }
        }

        if (!is.null(colnames(newdata_new))) {
            names_new <- colnames(newdata_new)
        } else {
            colnames(newdata_new) <- names_new
        }
    }

    while (any(names_new %in% colnames(x))) {
        names_new <- paste0(names_new, "_", names_new)
    }

    x_new <- x[1, ]
    for (i in 1:length(names_new)) {
        x_new[names_new[i]] <- 0
    }
    dists(x_new) <- dists_new
    sds(x_new) <- c(sds(x_new), sds_new)

    base_model <- attr(x_new, "base", exact = TRUE)
    fit_base <- attr(x_new, "fit_base_raw", exact = TRUE)

    if (base_model == "sep" && length(fit_base) == 2) {
        x_new <- add_base(x_new,
            fit_s = fit_base$fit_s,
            fit_t = fit_base$fit_t,
            sep = TRUE
        )
    } else {
        fit_base$dists_base <- dists_new$h
        x_new <- add_base.mcgf(x_new, fit_base = fit_base)
    }

    if (model == "all") {
        fit_lagr <- attr(x, "fit_lagr_raw", exact = TRUE)
        fit_lagr$dists_lagr <- dists_new
        x_new <- add_lagr.mcgf(x_new, fit_lagr = fit_lagr)
    }

    n_var_all <- n_var + n_var_new
    cov_mat <- ccov.mcgf(x_new, model = model)

    if (!no_newdata_new) {
        newdata_all <- cbind(newdata, newdata_new)

        cov_mat_res <- cov_par(cov_mat,
            horizon = horizon,
            n_var = n_var_all, joint = TRUE
        )
    } else {
        newdata_all <- newdata

        ind_rm <- n_var_all * horizon +
            unlist(lapply(
                0:(lag - 1),
                function(x) (n_var + 1):n_var_all + n_var_all * x
            ))
        cov_mat_res <- cov_par(cov_mat[-ind_rm, -ind_rm],
            horizon = horizon,
            n_var = n_var_all, joint = TRUE
        )
    }

    dat <- rbind(
        as.matrix(newdata_all),
        matrix(nrow = horizon - 1, ncol = ncol(newdata_all))
    )
    dat <- stats::embed(dat, n_block)
    pred <- tcrossprod(
        dat[, -c(1:(horizon * ncol(newdata_all)))],
        cov_mat_res$weights
    )
    locations_name <- colnames(newdata)

    if (is.null(locations_name)) {
        if (is.null(colnames(x))) {
            locations_name <- 1:n_var
        } else {
            locations_name <- colnames(x)
        }
    }

    Y_pred <- array(NA,
        dim = c(nrow(newdata_all), n_var_all, horizon),
        dimnames = list(
            rownames(newdata_all),
            c(locations_name, names_new),
            paste0("Horizon ", horizon:1)
        )
    )

    for (i in 1:horizon) {
        ind_y <- (n_block - i + 1):nrow(newdata_all)
        ind_pred <- 1:(nrow(dat) - horizon + i)
        Y_pred[ind_y, , i] <-
            pred[ind_pred, (1 + (i - 1) * n_var_all):(i * n_var_all)]
    }

    if (interval) {
        alpha <- (1 - level) / 2
        moe <- sqrt(diag(cov_mat_res$cov_curr)) *
            stats::qnorm(alpha, lower.tail = FALSE)

        lower <- upper <- array(NA,
            dim = dim(Y_pred),
            dimnames = dimnames(Y_pred)
        )
        for (i in 1:horizon) {
            moe_i <- moe[(1 + (i - 1) * n_var_all):(i * n_var_all)]
            lower[, , i] <- sweep(Y_pred[, , i], 2, moe_i)
            upper[, , i] <- sweep(Y_pred[, , i], 2, moe_i, "+")
        }

        Y_pred <- Y_pred[, , horizon:1]
        lower <- lower[, , horizon:1]
        upper <- upper[, , horizon:1]
        return(list(fit = Y_pred, lower = lower, upper = upper))
    } else {
        Y_pred <- Y_pred[, , horizon:1]
        return(Y_pred)
    }
}

#' Obtain kriging forecasts for new locations for an `mcgf_rs` object.
#'
#' @param x An `mcgf_rs` object.
#' @param newdata A data.frame with the same column names as `x`. If `newdata`
#' is missing the forecasts at the original data points are returned.
#' @param locations_new A matrix of data.frame of 2D points of new locations,
#' first column longitude, second column latitude, both in decimal degrees.
#' Supply only if `x` contains `locations`. Required when `dists_new_ls` is not
#' supplied.
#' @param dists_new_ls List of signed distance matrices (vectors) with names `h`,
#' `h1`, and 'h2' for all locations(and for each regime), with new locations
#' in the end. Each matrix must have the same number of columns. Required when
#' `locations_new` is not supplied.
#' @param newdata_new Optional; a data.frame with the same number of rows as
#' `newdata`. It contains the data of the new locations.
#' @param sds_new_ls List of the standard deviations of the new locations for
#' each regime. Format must be the same as the output from [`sds.mcgf_rs()`].
#' Default is 1 for all regimes.
#' @param newlabel A vector of new regime labels.
#' @param model Which model to use. One of `all`, `base`, or `empirical`.
#' @param interval Logical; if TRUE, prediction intervals are computed.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for
#' the intervals (if any) to be calculated. Used when `interval = TRUE`
#' @param soft Logical; if true, soft forecasts (and bounds) are produced.
#' @param prob Matrix with simplex rows. Number of columns must be the same as
#' unique labels in `x`.
#' @param dists_new_base Optional, list of distance matrices for the base
#' model. Used when the base model is non-regime switching. Default is `h` from
#' the first list of `dists_new_ls`.
#' @param ... Additional arguments.
#'
#' @return A list of kriging forecasts (and prediction intervals) for all
#' locations.
#' @export
#'
#' @details
#' It produces simple kriging forecasts for a zero-mean mcgf for new locations
#' given theri coordinates or relative distances. It supports kriging for the
#' `base` model and the `all` model which is the general stationary model with
#' the base and Lagrangian model from `x`.
#'
#' Users can either supply the coordinates via `locations_new`, or a list of
#' distance for all locations via `dists_new_ls`, with new locations at the
#' end. `dists_new_ls` will be used to calculate the new covariance matrices.
#' When `locations_new` is used, make sure `x` contains the attribute
#' `locations` of the coordinates of the old locations. When `dists_new_ls` is
#' used, it should be a list of a list of signed distance matrices of the same
#' dimension, where each row corresponds to the relative distances between a new
#' location and old locations in the same order as they appear in `x`. If only
#' one list is provided, it will be used for all regimes.
#'
#' When `soft = TRUE`, `prob` will be used to compute the soft forecasts
#' (weighted forecasts). The number of columns must match the number of unique
#' levels in `x`. The column order must be the same as the order of regimes as
#' in `levels(attr(x, "label", exact = TRUE))`. If not all regimes are seen in
#' `newlabel`, then only relevant columns in `prob` are used.
#'
#' When `interval = TRUE`, confidence interval for each forecasts and each
#' horizon is given. Note that it does not compute confidence regions.
#'
#' @examples
#' data(sim2)
#' sim2_mcgf <- mcgf_rs(sim2$data,
#'     locations = sim2$locations,
#'     label = sim2$label
#' )
#' sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
#' sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)
#'
#' # Fit a regime-switching separable model
#' fit_sep <- fit_base(
#'     sim2_mcgf,
#'     lag_ls = 5,
#'     model_ls = "sep",
#'     par_init_ls = list(list(
#'         c = 0.00005,
#'         gamma = 0.5,
#'         a = 0.5,
#'         alpha = 0.5
#'     )),
#'     par_fixed_ls = list(c(nugget = 0))
#' )
#'
#' # Store the fitted separable models to 'sim2_mcgf'
#' sim2_mcgf <- add_base(sim2_mcgf, fit_base_ls = fit_sep)
#'
#' # Calculate the simple kriging predictions and intervals for all locations
#' locations_new <- rbind(c(-110, 55), c(-109, 54))
#' sim2_krige <- krige_new(sim2_mcgf,
#'     locations_new = locations_new,
#'     model = "base", interval = TRUE
#' )
#' @family functions on fitting an mcgf_rs
krige_new.mcgf_rs <- function(x, newdata = NULL, locations_new = NULL,
                              dists_new_ls = NULL, newdata_new = NULL,
                              sds_new_ls = 1, newlabel,
                              soft = FALSE, prob, dists_new_base,
                              model = c("all", "base"),
                              interval = FALSE, level = 0.95, ...) {
    model <- match.arg(model)
    dots <- list(...)

    no_newdata <- ifelse(is.null(newdata), TRUE, FALSE)
    no_newdata_new <- ifelse(is.null(newdata_new), TRUE, FALSE)

    no_locations_new <- ifelse(is.null(locations_new), TRUE, FALSE)
    no_dists_new_ls <- ifelse(is.null(dists_new_ls), TRUE, FALSE)

    if (model == "base") {
        if (is.null(attr(x, "base", exact = T))) {
            stop("Base model missing from `x`.", call. = FALSE)
        }
    } else {
        if (is.null(attr(x, "lagr", exact = T))) {
            stop("Lagrangian model missing from `x`.", call. = FALSE)
        }
    }

    if (model == "base" && !attr(x, "base_rs", exact = TRUE)) {
        x_base <- x
        attr(x_base, "lag") <- attr(x, "lag")[[1]]
        attr(x_base, "sds") <- attr(x, "sds")$sds

        return(krige_new.mcgf(
            x = x_base, newdata = newdata,
            locations_new = locations_new,
            dists_new = dists_new_ls[[1]],
            newdata_new = newdata_new,
            sds_new = sds_new_ls[[1]], model = model, interval = interval,
            level = level, ...
        ))
    }

    if (model == "all" && !attr(x, "lagr_rs", exact = TRUE)) {
        x_lagr <- x
        attr(x_lagr, "lag") <- attr(x, "lag")[[1]]
        attr(x_lagr, "sds") <- attr(x, "sds")$sds

        return(krige.mcgf(
            x = x_lagr, newdata = newdata, locations_new = locations_new,
            dists_new = dists_new_ls[[1]], newdata_new = newdata_new,
            sds_new = sds_new_ls[[1]], model = model, interval = interval,
            level = level, ...
        ))
    }

    if (no_locations_new && no_dists_new_ls) {
        stop("must provide either `locations_new` or `dists_new_ls`.",
            call. = FALSE
        )
    }

    if (!no_locations_new && !no_dists_new_ls) {
        stop("do not provide both `locations_new` or `dists_new_ls`.",
            call. = FALSE
        )
    }

    if (!no_locations_new) {
        locations <- attr(x, "locations", exact = TRUE)

        if (min(locations_new[, 1]) < min(locations[, 1]) ||
            max(locations_new[, 1]) > max(locations[, 1])) {
            message("Longitude for new locations is outside of the grid.")
        }
        if (min(locations_new[, 2]) < min(locations[, 2]) ||
            max(locations_new[, 2]) > max(locations[, 2])) {
            message("Latitude for new locations is outside of the grid.")
        }

        if (is.null(locations)) {
            stop("`locations` not found in 'x', supply 'dists_new_ls' instead.",
                call. = FALSE
            )
        }

        if (is.vector(locations_new)) {
            if (length(locations_new) != 2) {
                stop("incorrect dimension for `locations_new`", call. = FALSE)
            }
            locations_new <- matrix(locations_new, nrow = 1)
        }
        n_var_new <- nrow(locations_new)
        names_new <- rownames(locations_new)

        if (!no_newdata_new) {
            if (ncol(newdata_new) != n_var_new) {
                stop("number of columns of `newdata_new` must be the same as ",
                    "the number of rows of `locations`",
                    call. = FALSE
                )
            }
        }

        longlat <- attr(x, "longlat", exact = TRUE)
        origin <- attr(x, "origin", exact = TRUE)
        dists_new_ls <- list(find_dists_new(locations, locations_new[, 1:2],
            longlat = longlat, origin = origin
        ))
    }

    lag_ls <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    n_var <- ncol(dists(x)$h)

    label <- attr(x, "label", exact = TRUE)
    lvs <- levels(label)
    n_regime <- length(lvs)

    if (!no_dists_new_ls) {
        if (!is.list(dists_new_ls)) {
            stop("`dists_new_ls` must be a list.", call. = FALSE)
        }
        if (!length(dists_new_ls) %in% c(1, n_regime)) {
            stop("length of `dists_new_ls` must be 1 or ", n_regime, ".",
                call. = FALSE
            )
        }

        dists_new_ls_dim <- lapply(
            dists_new_ls,
            function(x) lapply(x, dim)
        )
        if (length(table(unlist(dists_new_ls_dim))) != 1) {
            stop("dimensions must be the same for all elements in ",
                "`dists_new_ls`,",
                call. = FALSE
            )
        }

        if (any(!c("h", "h1", "h2") %in% names(dists_new_ls[[1]]))) {
            stop("elements in `dists_new_ls` must contain 'h', 'h1' and 'h2',",
                call. = FALSE
            )
        }
        if (any(dim(dists_new_ls[[1]]$h) != dim(dists_new_ls[[1]]$h1))) {
            stop("dimensions differ for 'h' and 'h1' in `dists_new_ls`.",
                call. = FALSE
            )
        }
        if (any(dim(dists_new_ls[[1]]$h1) != dim(dists_new_ls[[1]]$h2))) {
            stop("dimensions differ for 'h1' and 'h2' in `dists_new_ls`.",
                call. = FALSE
            )
        }
        if (is.array(dists_new_ls[[1]]$h1) &&
            !is.matrix(dists_new_ls[[1]]$h1)) {
            if (length(dim(dists_new_ls[[1]]$h1)) != 3) {
                stop("'h1' and 'h2' must be matrices or 3D arrays in ",
                    "`dists_new_ls`.",
                    call. = FALSE
                )
            }
        }

        n_var_new <- nrow(dists_new_ls[[1]]$h) - n_var
        names_new <-
            rownames(dists_new_ls[[1]]$h)[(n_var + 1):nrow(dists_new_ls[[1]]$h)]
    }

    if (length(dists_new_ls) == 1) {
        dists_new_ls <- rep(dists_new_ls, n_regime)
    }

    if (is.null(names_new)) {
        names_new <- paste0("New_", seq_len(n_var_new))
    }

    if (!is.list(sds_new_ls)) {
        if (length(sds_new_ls) == 1) {
            sds_new_ls <- rep(sds_new_ls, n_var_new)
        } else if (length(sds_new_ls) != n_var_new) {
            stop("`sds_new_ls` must be a list of length 1 or ",
                n_var_new,
                ".",
                call. = FALSE
            )
        }
        sds_new_ls <- rep(list(sds_new_ls), 2)
        sds_new_ls[[2]] <- rep(sds_new_ls[2], n_regime)
    } else if (is.list(sds_new_ls) && length(sds_new_ls) == 1) {
        if (length(sds_new_ls[[1]]) == 1) {
            sds_new_ls <- rep(sds_new_ls[[1]], n_var_new)
            sds_new_ls <- rep(list(sds_new_ls), 2)
            sds_new_ls[[2]] <- rep(sds_new_ls[2], n_regime)
        } else if (length(sds_new_ls[[1]]) != n_var_new) {
            stop("`sds_new_ls` must be a list of length 1 or ",
                n_var_new,
                ".",
                call. = FALSE
            )
        } else {
            sds_new_ls <- rep(sds_new_ls, 2)
            sds_new_ls[[2]] <- rep(sds_new_ls[2], n_regime)
        }
    } else {
        if (length(sds_new_ls) != 2) {
            stop("`sds_new_ls` must be a list of length 2.", call. = FALSE)
        } else if (length(sds_new_ls[[2]]) != n_regime) {
            stop("second element in `sds_new_ls` must be a list of length ",
                n_regime,
                ".",
                call. = FALSE
            )
        } else {
            if ((length(sds_new_ls[[1]]) != n_var_new) ||
                (length(sds_new_ls[[2]]) != n_regime)) {
                stop("incompatible `sds_new_ls` dimensions.", call. = FALSE)
            }
            for (i in 1:n_regime) {
                if (length(sds_new_ls[[2]][[i]]) != n_var_new) {
                    stop("incompatible `sds_new_ls` dimensions.",
                        call. = FALSE
                    )
                }
            }
        }
    }
    names(sds_new_ls) <- c("sds", "sds_rs")
    names(sds_new_ls[[1]]) <- names_new
    names(sds_new_ls[[2]]) <- paste("Regime", 1:n_regime)
    for (i in 1:n_regime) {
        names(sds_new_ls[[2]][[i]]) <- names_new
    }
    x_sds <- sds(x)
    sds_new_ls$sds <- c(x_sds$sds, sds_new_ls$sds)
    for (i in 1:n_regime) {
        sds_new_ls$sds_rs[[i]] <- c(
            x_sds$sds_rs[[i]],
            sds_new_ls$sds_rs[[i]]
        )
    }

    if (!no_newdata) {
        if (NCOL(newdata) != ncol(x)) {
            stop("unmatching number of columns for `newdata`.", call. = FALSE)
        }
        if (NROW(newdata) < max(unlist(lag_ls))) {
            stop("number of rows in `newdata` must be higher than `lag` ",
                max(unlist(lag_ls)), ".",
                call. = FALSE
            )
        }

        if (length(newlabel) != NROW(newdata)) {
            stop("lenght of `newlabel` must equal to `nrow(newdata)`.",
                call. = FALSE
            )
        }

        newlabel <- as.factor(newlabel)

        if (any(!(levels(newlabel) %in% lvs))) {
            stop("unknown levels in `newlabel.`", call. = FALSE)
        }

        label <- newlabel
    } else {
        newdata <- as.matrix(x)
    }

    if (!no_newdata_new) {
        if (ncol(newdata_new) != n_var_new) {
            stop("number of columns of `newdata_new` does not math with ",
                "`dists_new` or `locations_new`",
                call. = FALSE
            )
        }
        if (no_newdata) {
            if (nrow(newdata_new) != nrow(x)) {
                stop("number of rows do not match for `newdata_new` and `x`.",
                    call. = FALSE
                )
            }
        } else {
            if (nrow(newdata_new) != nrow(newdata)) {
                stop("number of rows do not match for `newdata_new` and ",
                    "`newdata`.",
                    call. = FALSE
                )
            }
        }

        if (!is.null(colnames(newdata_new))) {
            names_new <- colnames(newdata_new)
        } else {
            colnames(newdata_new) <- names_new
        }
    }

    while (any(names_new %in% colnames(x))) {
        names_new <- paste0(names_new, "_", names_new)
    }

    if (soft) {
        if (missing(prob)) {
            stop("must provide probabilities for soft forecasting.",
                call. = FALSE
            )
        }

        prob <- as.matrix(prob)

        if (ncol(prob) != length(lvs)) {
            stop("number of columns in `prob` must the same as the number of ",
                "unique levels in `x`.",
                call. = FALSE
            )
        }

        if (nrow(prob) != nrow(newdata)) {
            if (!no_newdata) {
                stop("number of rows in `prob` must be the same as that of ",
                    "`newdata`.",
                    call. = FALSE
                )
            } else {
                stop("number of rows in `prob` must be the same as that of ",
                    "`x`.",
                    call. = FALSE
                )
            }
        }

        new_lvs <- levels(label)

        prob[, which(!lvs %in% new_lvs)] <- 0

        if (ncol(prob) == 1) {
            soft <- FALSE
        } else {
            prob <- prob / rowSums(prob)
        }
    }

    x_new <- x[1, ]
    for (i in 1:length(names_new)) {
        x_new[names_new[i]] <- 0
    }
    dists(x_new) <- dists_new_ls[[1]]
    sds(x_new) <- sds_new_ls

    if (model == "base") {
        if (!attr(x, "base_rs", exact = TRUE)) {
            stop("`x` is not regime-switching, consider using an 'mcgf' object",
                call. = FALSE
            )
        }
    }

    if (!missing(dists_new_base)) {
        if (any(dim(dists_new_base) != dim(dists_new_ls[[1]]$h))) {
            stop("dimension of `dists_new_base` must be ",
                paste(dim(dists_new_ls[[1]]$h), collapse = " x "), ".",
                call. = FALSE
            )
        }
    } else {
        dists_new_base <- dists_new_ls[[1]]$h
    }

    fit_base <- attr(x, "fit_base_raw", exact = TRUE)
    if (attr(x, "base_rs", exact = TRUE)) {
        for (i in 1:n_regime) {
            fit_base[[i]]$dists_base <- dists_new_ls[[i]]$h
        }
        x_new <- add_base(x_new, fit_base_ls = fit_base)
    } else {
        fit_base$dists_base <- dists_new_base
        fit_base <- list(fit_base)
        fit_base$rs <- FALSE
        x_new <- add_base(x_new, fit_base_ls = fit_base)
    }

    if (model == "all") {
        if (!attr(x, "lagr_rs", exact = TRUE)) {
            stop("`x` is not regime-switching, consider using an 'mcgf' object",
                call. = FALSE
            )
        }

        fit_lagr <- attr(x, "fit_lagr_raw", exact = TRUE)
        for (i in 1:n_regime) {
            fit_lagr[[i]]$dists_lagr <- dists_new_ls[[i]]
        }
        x_new <- add_lagr(x_new, fit_lagr = fit_lagr)
    }

    n_var_all <- n_var + n_var_new
    cov_mat_ls <- ccov(x_new, model = model)

    if (!no_newdata_new) {
        newdata_all <- cbind(newdata, newdata_new)

        cov_mat_res <- lapply(cov_mat_ls, cov_par,
            horizon = horizon,
            n_var = n_var_all, joint = TRUE
        )
    } else {
        newdata_all <- newdata

        cov_mat_res <- vector("list", n_regime)
        for (i in 1:n_regime) {
            lag <- lag_ls[[i]]
            ind_rm <- n_var_all * horizon +
                unlist(lapply(
                    0:(lag - 1),
                    function(x) {
                        (n_var + 1):n_var_all + n_var_all * x
                    }
                ))
            cov_mat_res[[i]] <-
                cov_par(
                    cov_mat_ls[[i]][-ind_rm, -ind_rm],
                    horizon = horizon,
                    n_var = n_var_all,
                    joint = TRUE
                )
        }
    }

    locations_name <- colnames(x)

    if (is.null(locations_name)) {
        if (is.null(colnames(x))) {
            locations_name <- 1:n_var
        } else {
            locations_name <- colnames(x)
        }
    }

    Y_pred <- array(NA,
        dim = c(nrow(newdata_all), n_var_all, horizon),
        dimnames = list(
            rownames(newdata_all),
            c(locations_name, names_new),
            paste0("Horizon ", horizon:1)
        )
    )

    if (interval) {
        alpha <- (1 - level) / 2
        lower <- upper <- array(NA,
            dim = dim(Y_pred),
            dimnames = dimnames(Y_pred)
        )
        cv <- stats::qnorm(alpha, lower.tail = FALSE)
    }

    pred <- vector("list", n_regime)

    for (n in 1:n_regime) {
        lag_max <- lag_ls[[n]] + horizon - 1
        n_block <- lag_max + 1

        dat <- rbind(
            as.matrix(newdata_all),
            matrix(nrow = horizon - 1, ncol = ncol(newdata_all))
        )
        dat <- stats::embed(dat, n_block)
        pred[[n]] <- tcrossprod(
            dat[, -c(1:(horizon * ncol(newdata_all)))],
            cov_mat_res[[n]]$weights
        )
    }

    if (soft) {
        Y_pred_ls <- rep(list(Y_pred), n_regime)

        for (n in 1:n_regime) {
            for (i in 1:horizon) {
                ind_y <- (n_block - i + 1):nrow(newdata_all)
                ind_pred <- 1:(nrow(dat) - horizon + i)

                Y_pred_ls[[n]][ind_y, , i] <-
                    pred[[n]][ind_pred, (1 + (i - 1) * n_var_all):(i * n_var_all)] *
                        prob[ind_y, n]
            }
        }

        Y_pred <- Reduce("+", Y_pred_ls)

        if (interval) {
            sds_ls <- lapply(cov_mat_res, function(x) sqrt(diag(x$cov_curr)))
            sds_ls <- lapply(
                sds_ls,
                function(x) {
                    matrix(x,
                        nrow = nrow(prob),
                        ncol = n_var_all * horizon,
                        byrow = T
                    )
                }
            )
            sds_ls <- Map(
                function(x, i) x * prob[, i],
                sds_ls, seq_len(ncol(prob))
            )
            moe <- Reduce("+", sds_ls) * cv

            for (i in 1:horizon) {
                moe_i <- moe[, (1 + (i - 1) * n_var_all):(i * n_var_all)]
                lower[, , i] <- Y_pred[, , i] - moe_i
                upper[, , i] <- Y_pred[, , i] + moe_i
            }
        }
    } else {
        for (n in 1:n_regime) {
            for (i in 1:horizon) {
                ind_i <- (n_block - i + 1):nrow(newdata_all)
                ind_y <- ind_i[label[ind_i] == lvs[[n]]]

                ind_pred <- 1:(nrow(dat) - horizon + i)
                ind_pred <- ind_pred[label[ind_i] == lvs[[n]]]

                Y_pred[ind_y, , i] <-
                    pred[[n]][ind_pred, (1 + (i - 1) * n_var_all):(i * n_var_all)]

                if (interval) {
                    moe <- sqrt(diag(cov_mat_res[[n]]$cov_curr)) * cv
                    moe_i <- moe[(1 + (i - 1) * n_var_all):(i * n_var_all)]
                    lower[ind_y, , i] <-
                        sweep(Y_pred[ind_y, , i], 2, moe_i)
                    upper[ind_y, , i] <-
                        sweep(Y_pred[ind_y, , i], 2, moe_i, "+")
                }
            }
        }
    }

    Y_pred <- Y_pred[, , horizon:1]

    if (interval) {
        lower <- lower[, , horizon:1]
        upper <- upper[, , horizon:1]

        return(list(fit = Y_pred, lower = lower, upper = upper))
    } else {
        return(Y_pred)
    }
}
