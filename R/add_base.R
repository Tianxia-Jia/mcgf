#' Generic function for adding a base model
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @details
#' Refer to [`add_base.mcgf()`] for more details.
#'
#' @export
add_base <- function(x, ...) {
    UseMethod("add_base")
}

#' Add base model outputted from [`fit_base()`] to an `mcgf` or `mcgf_rs`
#' object.
#'
#' @name add_base.mcgf
#'
#' @param x An `mcgf` or `mcgf_rs` object.
#' @param fit_base Output from the [`fit_base()`] function.
#' @param fit_s Pure spatial model outputted from the [`fit_base()`] function.
#' Used only when `sep = TRUE`.
#' @param fit_t Pure temporal model outputted from the [`fit_base()`] function.
#' Used only when `sep = TRUE`.
#' @param sep Logical; TRUE if spatial and temporal models are fitted
#' separately.
#' @param old Logical; TRUE if the old base model needs to be kept.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the base model.
#' @export
#'
#' @details
#' After fitting the base model by [`fit_base()`], the results can be added to
#' `x` by [`add_base()`]. To supply the base model directly, use [`base<-`] to
#' add the base model; the value must contain `model`, `par_base`, `cor_base`,
#' `lag`, and `horizon`.
#'
#' @family {functions related to model fitting}
add_base.mcgf <-
    function(x,
             fit_base,
             fit_s,
             fit_t,
             sep = F,
             old = FALSE,
             ...) {

    if (old) {

        attr(x, "base_old") <- attr(x, "base", exact = TRUE)
        attr(x, "base_res_old") <- attr(x, "base_res", exact = TRUE)
        attr(x, "horizon_old") <- attr(x, "horizon", exact = TRUE)

        if (sep) {
            lag_new <- fit_t$lag
            horizon_new <- fit_t$horizon

        } else {
            lag_new <- fit_base$lag
            horizon_new <- fit_base$horizon
        }

        if (!is.null(attr(x, "lag", exact = TRUE)) &&
            attr(x, "lag", exact = TRUE) != lag_new)
            warning("unmatching `lag` for old and new base models.")

        if (!is.null(attr(x, "horizon", exact = TRUE)) &&
            attr(x, "horizon", exact = TRUE) != horizon_new)
            warning("unmatching `horizon` for old and new base models.")
    }

    if (sep) {

        if (missing(fit_s) || missing(fit_t))
            stop("must give `fit_s` and `fit_t`.")

        par_s <- as.list(fit_s$fit$par)
        names(par_s) <- fit_s$par_names
        par_s <- c(par_s, fit_t$par_fixed)

        par_t <- as.list(fit_t$fit$par)
        names(par_t) <- fit_t$par_names
        par_t <- c(par_t, fit_t$par_fixed)

        par_base <- c(par_s, par_t)

        base_fn <- "..cor_sep"
        lag_max <- fit_t$lag + fit_t$horizon - 1

    } else {

        if (!(fit_base$model %in% c("sep", "fs")))
            stop('base model must be `sep` or `fs`.')

        par_base <- as.list(fit_base$fit$par)
        names(par_base) <- fit_base$par_names
        par_base <- c(par_base, fit_base$par_fixed)

        base_fn <- switch(fit_base$model,
                          sep = "..cor_sep",
                          fs = ".cor_fs")

        lag_max <- fit_base$lag + fit_base$horizon - 1
    }

    h_u_ar <- to_ar(h = dists(x)$h, lag_max = lag_max)
    par_base_other <- list(h = h_u_ar$h_ar, u = h_u_ar$u_ar)

    cor_base <- do.call(base_fn, c(par_base, par_base_other))

    if (sep) {

        base_res <- list(
            par_base = par_base,
            fit_base = list(spatial = fit_s$fit, temporal = fit_t$fit),
            method_base = c(spatial = fit_s$method, temporal = fit_t$method),
            optim_fn = c(spatial = fit_s$optim_fn, temporal = fit_t$optim_fn),
            cor_base = cor_base,
            par_fixed = c(fit_s$par_fixed, fit_t$par_fixed),
            dots = list(spatial = fit_s$dots, temporal = fit_t$dots)
        )
        attr(x, "base") <- "sep"
        attr(x, "base_res") <- base_res
        attr(x, "lag") <- fit_t$lag
        attr(x, "horizon") <- fit_t$horizon

    } else {

        base_res <- list(
            par_base = par_base,
            fit_base = fit_base$fit,
            method_base = fit_base$method,
            optim_fn = fit_base$optim_fn,
            cor_base = cor_base,
            par_fixed = fit_base$par_fixed,
            dots = fit_base$dots
        )
        dists_base <- fit_base$dists_base
        base_res <- c(base_res, dists_base = dists_base)

        attr(x, "base") <- fit_base$model
        attr(x, "base_res") <- base_res
        attr(x, "lag") <- fit_base$lag
        attr(x, "horizon") <- fit_base$horizon
    }
    return(x)
}

#' @rdname add_base.mcgf
#'
#' @param value A list containing the base model as well as its parameters. It
#' must contains `model`, `par_base`, `cor_base`, `lag`, and `horizon`.
#' @export
`base<-` <- function(x, value) {

    if (any(! c("model", "par_base", "cor_base", "lag", "horizon") %in%
            names(value)))
        stop("`value` must contain `model`, `par_base`, `cor_base`, `lag`",
             ", and `horizon`.")

    if(is.null(attr(x, "base", exact = TRUE)))
        message("Overwriting the existing base model.")

    base_res <- list(
        par_base = value$par_base,
        fit_base = value$fit,
        method_base = value$method,
        optim_fn = value$optim_fn,
        cor_base = value$cor_base,
        par_fixed = value$par_fixed,
        dots = value$dots
    )
    dists_base <- value$dists_base
    base_res <- c(base_res, dists_base = dists_base)

    attr(x, "base") <- value$model
    attr(x, "base_res") <- base_res
    attr(x, "lag") <- value$lag
    attr(x, "horizon") <- value$horizon
    return(x)
}
