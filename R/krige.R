#' Generic function for displaying fitted models for mcgf objects
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


#' Obtain kriging forecasts for an mcgf object.
#'
#' @param x An mcgf object.
#' @param newdata A data.frame with the same column names as `x`. If `newdata`
#' is missing the forecasted values at the original data points are returned.
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
            if (is.null(lag))
                stop("please provide `lag` for the empirical model.")
        }
        if (is.null(horizon)) {
            horizon <- dots$horizon
            if (is.null(horizon))
                stop("please provide `horizon` for the empirical model.")
        }

    } else {

        if (model == "base") {
            if (is.null(attr(x, "base", exact = T)))
                stop("Base model missing from `x`.")

        } else {
            if (is.null(attr(x, "lagr", exact = T)))
                stop("Lagrangian model missing from `x`.")
        }
    }

    lag_max <- lag + horizon - 1
    n_blcok <- lag_max + 1
    n_var <- ncol(dists(x)$h)
    cov_mat <- cov_par(ccov(x, model = model), horizon = horizon, n_var = n_var,
                       joint = TRUE)

    if (!missing(newdata)) {

        if (NCOL(newdata) != ncol(x))
            stop("unmatching number of columns for `newdata`.")
        if (NROW(newdata) < lag)
            stop("number of rows in `newdata` must be higher than `lag` ",
                 lag, ".")

        dat <- stats::embed(as.matrix(newdata), n_blcok)
    } else {
        dat <- stats::embed(as.matrix(x), n_blcok)
    }

    #Y <- dat[, 1:(horizon * n_var)]
    pred <- dat[, -c(1:(horizon * n_var))] %*% t(cov_mat$weights)
    n_pred <- nrow(pred)

    Y_pred <- array(NA, dim = c(dim(x), horizon),
                    dimnames = list(rownames(x),
                                    colnames(x),
                                    paste0("Horizon ", horizon:1)))

    for (i in 1:horizon) {
        Y_pred[(n_blcok - i + 1):(nrow(x) - i + 1), , i] <-
            pred[, (1 + (i - 1) * n_var):(i * n_var)]
    }

    if (interval) {

        alpha <- (1 - level) / 2
        moe <- sqrt(diag(cov_mat$cov_curr)) *
            stats::qnorm(alpha, lower.tail = FALSE)

        lower <- upper <- array(NA, dim = dim(Y_pred),
                                dimnames = dimnames(Y_pred))
        for (i in 1:horizon) {
            moe_hori <- moe[(1 + (i - 1) * n_var):(i * n_var)]
            lower[, , i] <- sweep(Y_pred[, , i], 2, moe_hori)
            upper[, , i] <- sweep(Y_pred[, , i], 2, moe_hori, "+")
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

