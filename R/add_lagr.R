#' Generic function for adding a Lagrangian model
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return lagr model and its parameter
#' @export
#' @family {functions related to model fitting}
add_lagr <- function(x, ...) {
    UseMethod("add_lagr")
}

#' Add lagr model outputted from [`fit_lagr()`] to a mcgf object.
#'
#' @name add_lagr.mcgf
#'
#' @param x An mcgf object.
#' @param fit_lagr Output from the [`fit_lagr()`] function.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the Lagrangian model.
#' @export
#' @family {functions related to model fitting}
add_lagr.mcgf <- function(x, fit_lagr, ...) {

    par_lagr <- as.list(fit_lagr$fit$par)
    names(par_lagr) <- fit_lagr$par_names
    par_lagr <- c(par_lagr, fit_lagr$par_fixed)

    lagrangian <- fit_lagr$model
    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_max <- lag + horizon - 1

    cor_base <- attr(x, "base_res", exact = TRUE)$cor_base
    u_ar <- to_ar(h = dists(x)$h, lag_max = lag_max)$u_ar
    h1_ar <- to_ar(h = dists(x)$h1, lag_max = lag_max, u = FALSE)
    h2_ar <- to_ar(h = dists(x)$h2, lag_max = lag_max, u = FALSE)

    par_lagr_other <- list(
        cor_base = cor_base,
        lagrangian = lagrangian,
        h1 = h1_ar,
        h2 = h2_ar,
        u = u_ar
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

    attr(x, "lagr") <- fit_lagr$model
    attr(x, "lagr_res") <- lagr_res

    return(x)
}

#' @rdname add_lagr.mcgf
#'
#' @param value A list containing the lagr model as well as its parameters. It
#' must contains `model`, `par_lagr`, and `cor_lagr`.
#' @export
`lagr<-` <- function(x, value) {

    if (any(! c("model", "par_lagr", "cor_lagr") %in%
            names(value)))
        stop("`value` must contain `model`, `par_lagr`, `cor_lagr`.")

    if(is.null(attr(x, "lagr", exact = TRUE)))
        message("Overwriting the existing lagr model.")

    lagr_res <- list(
        par_lagr = value$par_lagr,
        fit_lagr = value$fit,
        method_lagr = value$method,
        optim_fn = value$optim_fn,
        cor_lagr = value$cor_lagr,
        par_fixed = value$par_fixed,
        dots = value$dots
    )

    attr(x, "lagr") <- value$model
    attr(x, "lagr_res") <- lagr_res
    return(x)
}
