#' Generic function for adding a Lagrangian model
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return `x` with the newly added Lagrangian model.
#' @export
#'
#' @details
#' Refer to [`add_lagr.mcgf()`] and [`add_lagr.mcgf_rs()`] for more details.
add_lagr <- function(x, ...) {
    UseMethod("add_lagr")
}

#' Add lagr model outputted from [`fit_lagr()`] to a `mcgf` object.
#'
#' @name add_lagr.mcgf
#'
#' @param x An `mcgf` object.
#' @param fit_lagr Output from the [`fit_lagr()`] function.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the Lagrangian model.
#' @export
#'
#' @examples
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
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
#' model(sim1_mcgf, old = TRUE)
#' @family functions related to fitting an `mcgf` object
add_lagr.mcgf <- function(x, fit_lagr, ...) {
    par_lagr <- as.list(fit_lagr$fit$par)
    names(par_lagr) <- fit_lagr$par_names
    par_lagr <- c(par_lagr, fit_lagr$par_fixed)

    lagrangian <- fit_lagr$model
    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_max <- lag + horizon - 1
    scale_time <- attr(x, "scale_time", exact = TRUE)

    if (!is.null(fit_lagr$dists_lagr)) {
        lagr_h <- fit_lagr$dists_base$h
        lagr_h1 <- fit_lagr$dists_base$h1
        lagr_h2 <- fit_lagr$dists_base$h2
    } else {
        lagr_h <- dists(x)$h
        lagr_h1 <- dists(x)$h1
        lagr_h2 <- dists(x)$h2
    }

    if (!is.matrix(lagr_h)) {
        lagr_h <- lagr_h[, , 1:(lag_max + 1)]
        lagr_h1 <- lagr_h1[, , 1:(lag_max + 1)]
        lagr_h2 <- lagr_h2[, , 1:(lag_max + 1)]
    }

    cor_base <- attr(x, "base_res", exact = TRUE)$cor_base
    u_ar <- to_ar(h = lagr_h, lag_max = lag_max)$u_ar
    h1_ar <- to_ar(h = lagr_h1, lag_max = lag_max, u = FALSE)
    h2_ar <- to_ar(h = lagr_h2, lag_max = lag_max, u = FALSE)

    par_lagr_other <- list(
        cor_base = cor_base,
        lagrangian = lagrangian,
        h1 = h1_ar,
        h2 = h2_ar,
        u = u_ar / scale_time
    )

    cor_lagr <- do.call("..cor_stat", c(par_lagr, par_lagr_other))

    lagr_res <- list(
        par_lagr = par_lagr,
        fit_lagr = fit_lagr$fit,
        method_lagr = fit_lagr$method,
        optim_fn = fit_lagr$optim_fn,
        cor_lagr = cor_lagr,
        par_fixed = fit_lagr$par_fixed,
        dots = fit_lagr$dots
    )

    dists_lagr <- fit_lagr$dists_lagr
    lagr_res <- c(lagr_res, list(dists_lagr = dists_lagr))

    attr(x, "lagr") <- fit_lagr$model
    attr(x, "lagr_res") <- lagr_res

    return(x)
}

#' Add lagr model outputted from [`fit_lagr()`] to a `mcgf_rs` object.
#'
#' @param x An `mcgf_rs` object.
#' @param fit_lagr_ls Output from the [`fit_lagr()`] function.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the Lagrangian model.
#' @export
#'
#' @details
# ‘ This function is equivalent to [`add_lagr.mcgf()`] for `mcgf_rs` objects.
# ‘
#' After fitting the Lagrangian model by [`fit_lagr()`], the results can be
#' added to `x` by [`add_base()`]. To supply the Lagrangian model directly,
#' use [`lagr<-`] to add the Lagrangian model; the value must contain the same
#' output as [`add_lagr.mcgf()`] or [`add_lagr.mcgf_rs()`].
#'
#' @examples
#' data(sim3)
#' sim3_mcgf <- mcgf_rs(sim3$data, dists = sim3$dists, label = sim3$label)
#' sim3_mcgf <- add_acfs(sim3_mcgf, lag_max = 5)
#' sim3_mcgf <- add_ccfs(sim3_mcgf, lag_max = 5)
#'
#' # Fit a fully symmetric model with known variables
#' fit_fs <- fit_base(
#'     sim3_mcgf,
#'     lag_ls = 5,
#'     model_ls = "fs",
#'     rs = FALSE,
#'     par_init_ls = list(list(beta = 0)),
#'     par_fixed_ls = list(list(
#'         nugget = 0,
#'         c = 0.05,
#'         gamma = 0.5,
#'         a = 0.5,
#'         alpha = 0.2
#'     ))
#' )
#'
#' # Set beta to 0 to fit a separable model with known variables
#' fit_fs[[1]]$fit$par <- 0
#'
#' # Store the fitted separable model to 'sim3_mcgf'
#' sim3_mcgf <- add_base(sim3_mcgf, fit_base_ls = fit_fs)
#'
#' # Fit a regime-switching Lagrangian model.
#' fit_lagr_rs <- fit_lagr(
#'     sim3_mcgf,
#'     model_ls = list("lagr_tri"),
#'     par_init_ls = list(
#'         list(v1 = -50, v2 = 50),
#'         list(v1 = 100, v2 = 100)
#'     ),
#'     par_fixed_ls = list(list(lambda = 0.2, k = 2))
#' )
#'
#' # Store the fitted Lagrangian model to 'sim3_mcgf'
#' sim3_mcgf <- add_lagr(sim3_mcgf, fit_lagr_ls = fit_lagr_rs)
#' model(sim3_mcgf)
#' @family functions related to fitting an `mcgf_rs` object
add_lagr.mcgf_rs <- function(x, fit_lagr_ls, ...) {
    if (!fit_lagr_ls$rs) {
        attr(x, "lag") <- attr(x, "lag")[[1]]
        x <- add_lagr.mcgf(x = x, fit_lagr = fit_lagr_ls[[1]], ...)
        attr(x, "lagr_rs") <- FALSE
        return(x)
    }

    lvs <- levels(attr(x, "label", exact = TRUE))
    n_regime <- length(lvs)

    lag_ls <- attr(x, "lag", exact = TRUE)
    if (length(lag_ls) == 1) {
        lag_ls <- rep(lag_ls, n_regime)
        names(lag_ls) <- paste0("Regime ", lvs)
        attr(x, "lag") <- lag_ls
    }

    lagr_res_ls <- lagr_model_ls <- vector("list", n_regime)
    names(lagr_res_ls) <- names(lagr_model_ls) <- paste0("Regime ", lvs)
    scale_time <- attr(x, "scale_time", exact = TRUE)

    for (i in 1:n_regime) {
        fit_lagr <- fit_lagr_ls[[i]]

        par_lagr <- as.list(fit_lagr$fit$par)
        names(par_lagr) <- fit_lagr$par_names
        par_lagr <- c(par_lagr, fit_lagr$par_fixed)

        lagrangian <- fit_lagr$model
        lag <- lag_ls[[i]]
        horizon <- attr(x, "horizon", exact = TRUE)
        lag_max <- lag + horizon - 1

        if (!is.null(fit_lagr$dists_lagr)) {
            lagr_h <- fit_lagr$dists_base$h
            lagr_h1 <- fit_lagr$dists_base$h1
            lagr_h2 <- fit_lagr$dists_base$h2
        } else {
            lagr_h <- dists(x)$h
            lagr_h1 <- dists(x)$h1
            lagr_h2 <- dists(x)$h2
        }

        if (!is.matrix(lagr_h)) {
            lagr_h <- lagr_h[, , 1:(lag_max + 1)]
            lagr_h1 <- lagr_h1[, , 1:(lag_max + 1)]
            lagr_h2 <- lagr_h2[, , 1:(lag_max + 1)]
        }

        if (any(attr(x, "base_rs", exact = TRUE))) {
            cor_base <- attr(x, "base_res", exact = TRUE)[[i]]$cor_base
        } else {
            cor_base <- attr(x, "base_res", exact = TRUE)$cor_base
        }

        u_ar <- to_ar(h = lagr_h, lag_max = lag_max)$u_ar
        h1_ar <- to_ar(h = lagr_h1, lag_max = lag_max, u = FALSE)
        h2_ar <- to_ar(h = lagr_h2, lag_max = lag_max, u = FALSE)

        par_lagr_other <- list(
            cor_base = cor_base,
            lagrangian = lagrangian,
            h1 = h1_ar,
            h2 = h2_ar,
            u = u_ar / scale_time
        )
        cor_lagr <- do.call("..cor_stat", c(par_lagr, par_lagr_other))

        lagr_res <- list(
            par_lagr = par_lagr,
            fit_lagr = fit_lagr$fit,
            method_lagr = fit_lagr$method,
            optim_fn = fit_lagr$optim_fn,
            cor_lagr = cor_lagr,
            par_fixed = fit_lagr$par_fixed,
            dots = fit_lagr$dots
        )
        dists_lagr <- fit_lagr$dists_lagr
        lagr_res <- c(lagr_res, list(dists_lagr = dists_lagr))

        lagr_res_ls[[i]] <- lagr_res
        lagr_model_ls[[i]] <- lagrangian
    }

    attr(x, "lagr") <- lagr_model_ls
    attr(x, "lagr_res") <- lagr_res_ls
    attr(x, "lagr_rs") <- fit_lagr_ls$rs

    return(x)
}

#' @rdname add_lagr.mcgf
#'
#' @param value A list containing the lagr model as well as its parameters. It
#' must contains `model`, `par_lagr`, and `cor_lagr`.
#' @export
`lagr<-` <- function(x, value) {
    if (any(!c("model", "par_lagr", "cor_lagr") %in%
        names(value))) {
        stop("`value` must contain `model`, `par_lagr`, `cor_lagr`.")
    }

    if (is.null(attr(x, "lagr", exact = TRUE))) {
        message("Overwriting the existing lagr model.")
    }

    lagr_res <- list(
        par_lagr = value$par_lagr,
        fit_lagr = value$fit,
        method_lagr = value$method,
        optim_fn = value$optim_fn,
        cor_lagr = value$cor_lagr,
        par_fixed = value$par_fixed,
        dots = value$dots
    )
    dists_lagr <- value$dists_lagr
    lagr_res <- c(lagr_res, list(dists_lagr = dists_lagr))

    attr(x, "lagr") <- value$model
    attr(x, "lagr_res") <- lagr_res
    attr(x, "lagr_rs") <- value$rs
    return(x)
}
