#' Generic function for computing kriging forecasts
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Kriging results of `x`
#' @export
#' @family {functions related to prediction}
krige <- function(x, ...) {
    UseMethod("krige")
}


#' Obtain kriging forecasts for an `mcgf` object.
#'
#' @param x An `mcgf` object.
#' @param newdata A data.frame with the same column names as `x`. If `newdata`
#' is missing the forecasts at the original data points are returned.
#' @param model Which model to use. One of `all`, `base`, and `empirical`.
#' @param interval Logical; if TRUE, prediction intervals are computed.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for
#' the intervals (if any) to be calculated. Used when `interval = TRUE`
#' @param ... Additional arguments. Give `lag` and `horizon` if they are not
#' defined in `x` for the `empirical` model.
#'
#' @return A list of kriging forecasts (and prediction intervals).
#' @export
#'
#' @details
#' It produces simple kriging forecasts for a zero-mean mcgf. It supports
#' kriging for the `empirical` model, the `base` model, and the `all` model
#' which is the general stationary model with the base and Lagrangian model
#' from `x`.
#'
#' When `interval = TRUE`, confidence interval for each forecasts and each
#' horizon is given. Note that it does not compute confidence regions.
#'
#' @family {functions related to prediction}
krige.mcgf <- function(x, newdata, model = c("all", "base", "empirical"),
                       interval = FALSE, level = 0.95, ...) {
    model <- match.arg(model)
    dots <- list(...)

    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)

    if (model == "empirical") {
        if (is.null(lag)) {
            lag <- dots$lag
            if (is.null(lag)) {
                stop("please provide `lag` for the empirical model.",
                    call. = FALSE
                )
            }
            attr(x, "lag") <- lag
        }
        if (is.null(horizon)) {
            horizon <- dots$horizon
            if (is.null(horizon)) {
                stop("please provide `horizon` for the empirical model.",
                    call. = FALSE
                )
            }
            attr(x, "horizon") <- horizon
        }
    } else if (model == "base") {
        if (is.null(attr(x, "base", exact = T))) {
            stop("Base model missing from `x`.", call. = FALSE)
        }
    } else {
        if (is.null(attr(x, "lagr", exact = T))) {
            stop("Lagrangian model missing from `x`.", call. = FALSE)
        }
    }

    lag_max <- lag + horizon - 1
    n_block <- lag_max + 1
    n_var <- ncol(dists(x)$h)
    cov_mat <- ccov.mcgf(x, model = model)
    cov_mat_res <- cov_par(cov_mat,
        horizon = horizon,
        n_var = n_var, joint = TRUE
    )

    if (!missing(newdata)) {
        if (NCOL(newdata) != ncol(x)) {
            stop("unmatching number of columns for `newdata`.", call. = FALSE)
        }
        if (NROW(newdata) < lag) {
            stop("number of rows in `newdata` must be higher than `lag` ",
                lag, ".",
                call. = FALSE
            )
        }

        x <- newdata
    }

    dat <- rbind(as.matrix(x), matrix(nrow = horizon - 1, ncol = ncol(x)))
    dat <- stats::embed(dat, n_block)
    pred <- dat[, -c(1:(horizon * n_var))] %*% t(cov_mat_res$weights)
    Y_pred <- array(NA,
        dim = c(dim(x), horizon),
        dimnames = list(
            rownames(x),
            colnames(x),
            paste0("Horizon ", horizon:1)
        )
    )

    for (i in 1:horizon) {
        ind_y <- (n_block - i + 1):nrow(x)
        ind_pred <- 1:(nrow(dat) - horizon + i)
        Y_pred[ind_y, , i] <- pred[ind_pred, (1 + (i - 1) * n_var):(i * n_var)]
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
            moe_i <- moe[(1 + (i - 1) * n_var):(i * n_var)]
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

#' Obtain kriging forecasts for an `mcgf_rs` object.
#'
#' @param x An `mcgf_rs` object.
#' @param newdata A data.frame with the same column names as `x`. If `newdata`
#' is missing the forecasts at the original data points are returned.
#' @param newlabel A vector of new regime labels.
#' @param model Which model to use. One of `all`, `base`, and `empirical`.
#' @param interval Logical; if TRUE, prediction intervals are computed.
#' @param level A numeric scalar between 0 and 1 giving the confidence level for
#' the intervals (if any) to be calculated. Used when `interval = TRUE`
#' @param soft Logical; if true, soft forecasts (and bounds) are produced.
#' @param prob Matrix with simplex rows. Number of columns must be the same as
#' unique labels in `x`.
#' @param ... Additional arguments. Give `lag` and `horizon` if they are not
#' defined in `x` for the `empirical` model.
#'
#' @return A list of kriging forecasts (and prediction intervals).
#' @export
#'
#' @details
#' It produces simple kriging forecasts for a zero-mean mcgf. It supports
#' kriging for the `empirical` model, the `base` model, and the `all` model
#' which is the general stationary model with the base and Lagrangian model
#' from `x`.
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
#' @family {functions related to prediction}
krige.mcgf_rs <- function(x, newdata, newlabel,
                          soft = FALSE,
                          prob,
                          model = c("all", "base", "empirical"),
                          interval = FALSE, level = 0.95, ...) {
    model <- match.arg(model)

    if (model == "base" & !attr(x, "base_rs", exact = TRUE)) {
        x_base <- x
        attr(x_base, "lag") <- attr(x, "lag")[[1]]
        attr(x_base, "sds") <- attr(x, "sds")$sds

        return(krige.mcgf(
            x = x_base, newdata = newdata, model = model,
            interval = interval, level = level, ...
        ))
    }

    if (model == "all" & !attr(x, "lagr_rs", exact = TRUE)) {
        x_lagr <- x
        attr(x_lagr, "lag") <- attr(x, "lag")[[1]]
        attr(x_lagr, "sds") <- attr(x, "sds")$sds

        return(krige.mcgf(
            x = x_lagr, newdata = newdata, model = model,
            interval = interval, level = level, ...
        ))
    }

    dots <- list(...)

    lag_ls <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    n_var <- ncol(dists(x)$h)

    label <- attr(x, "label", exact = TRUE)
    lvs <- levels(label)
    n_regime <- length(lvs)

    if (model == "empirical") {
        if (is.null(lag_ls)) {
            lag_ls <- dots$lag_ls

            if (is.null(lag_ls)) {
                stop("please provide `lag_ls` for the empirical model.",
                    call. = FALSE
                )
            }

            if (length(lag_ls) != n_regime) {
                stop("length of `lag_ls` must be ", n_regime, ".",
                    call. = FALSE
                )
            }
        }

        if (is.null(horizon)) {
            horizon <- dots$horizon
            if (is.null(horizon)) {
                stop("please provide `horizon` for the empirical model.",
                    call. = FALSE
                )
            }
        }
    } else if (model == "base") {
        if (is.null(attr(x, "base", exact = T))) {
            stop("Base model missing from `x`.")
        }
    } else {
        if (is.null(attr(x, "lagr", exact = T))) {
            stop("Lagrangian model missing from `x`.")
        }
    }

    cov_mat_ls <- ccov(x, model = model)
    cov_mat_res <- lapply(cov_mat_ls, cov_par,
        horizon = horizon,
        n_var = n_var, joint = TRUE
    )

    if (!missing(newdata)) {
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

        x <- newdata
        label <- newlabel
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

        if (nrow(prob) != nrow(x)) {
            if (!missing(newdata)) {
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

        if (nrow(prob) != nrow(x)) {
            stop("number of rows in `prob` must be the same as that of ",
                "`x`.",
                call. = FALSE
            )
        }

        new_lvs <- levels(label)

        prob[, which(!lvs %in% new_lvs)] <- 0

        if (ncol(prob) == 1) {
            soft <- FALSE
        } else {
            prob <- prob / rowSums(prob)
        }
    }

    Y_pred <- array(NA,
        dim = c(dim(x), horizon),
        dimnames = list(
            rownames(x),
            colnames(x),
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
            as.matrix(x),
            matrix(nrow = horizon - 1, ncol = ncol(x))
        )
        dat <- stats::embed(as.matrix(dat), n_block)

        pred[[n]] <- dat[, -c(1:(horizon * n_var))] %*%
            t(cov_mat_res[[n]]$weights)
    }

    if (soft) {
        Y_pred_ls <- rep(list(Y_pred), n_regime)

        for (n in 1:n_regime) {
            for (i in 1:horizon) {
                ind_y <- (n_block - i + 1):nrow(x)
                ind_pred <- 1:(nrow(dat) - horizon + i)

                Y_pred_ls[[n]][ind_y, , i] <-
                    pred[[n]][ind_pred, (1 + (i - 1) * n_var):(i * n_var)] *
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
                        ncol = n_var * horizon,
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
                moe_i <- moe[, (1 + (i - 1) * n_var):(i * n_var)]
                lower[, , i] <- Y_pred[, , i] - moe_i
                upper[, , i] <- Y_pred[, , i] + moe_i
            }
        }
    } else {
        for (n in 1:n_regime) {
            for (i in 1:horizon) {
                ind_i <- (n_block - i + 1):nrow(x)
                ind_y <- ind_i[label[ind_i] == lvs[[n]]]

                ind_pred <- 1:(nrow(dat) - horizon + i)
                ind_pred <- ind_pred[label[ind_i] == lvs[[n]]]

                Y_pred[ind_y, , i] <-
                    pred[[n]][ind_pred, (1 + (i - 1) * n_var):(i * n_var)]

                if (interval) {
                    moe <- sqrt(diag(cov_mat_res[[n]]$cov_curr)) * cv
                    moe_i <- moe[(1 + (i - 1) * n_var):(i * n_var)]
                    lower[ind_y, , i] <-
                        sweep(Y_pred[ind_y, , i], 2, moe_i)
                    upper[ind_y, , i] <-
                        sweep(Y_pred[ind_y, , i], 2, moe_i, "+")
                }
            }
        }
    }

    if (interval) {
        Y_pred <- Y_pred[, , horizon:1]
        lower <- lower[, , horizon:1]
        upper <- upper[, , horizon:1]

        return(list(fit = Y_pred, lower = lower, upper = upper))
    } else {
        Y_pred <- Y_pred[, , horizon:1]
        return(Y_pred)
    }
}
