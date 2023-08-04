#' Generic function for adding a base model
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Base model and its parameter
#' @export
#' @family {functions related to model fitting}
add_base <- function(x, ...) {
    UseMethod("add_base")
}

#' Add base model outputted from `fit_base` to a mcgf object.
#'
#' @name add_base.mcgf
#'
#' @param x An mcgf object.
#' @param fit_base Output from the `fit_base` function.
#' @param ... Additional arguments. Not in use.
#'
#' @return `x` with newly added attributes of the base model.
#' @export
#' @family {functions related to model fitting}
add_base.mcgf <- function(x, fit_base, ...) {

    if (!(fit_base$base %in% c("sep", "fs")))
        stop('base model must be `sep` or `fs`.')

    par_base <- list(fit_base$par)
    names(par_base) <- fit_base$par_names
    par_base <- c(par_base, fit_base$par_fixed)

    base_fn <- switch(fit_base$base,
                       sep = "..cor_sep",
                       fs = ".cor_fs")

    lag_max <- fit_base$lag + fit_base$horizon - 1

    h_u_ar <- to_ar(h = dists(x)$h, lag_max = lag_max)
    par_base_other <- list(h = h_u_ar$h_ar, u = h_u_ar$u_ar)

    cov_base <- do.call(base_fn, c(par_base, par_base_other))
    cov_base <- cor2cov_ar(cov_base, sds(x))
    cov_base <- cov_joint(cov_base)

    base_res <- list(
        par_base = par_base,
        fit_base = fit_base$fit,
        cov_base = cov_base,
        lag_base = fit_base$lag
    )

    attr(x, "base") <- fit_base$model
    attr(x, "base_res") <- base_res
    attr(x, "horizon") <- fit_base$horizon
    return(x)
}

#' @rdname add_base.mcgf
#'
#' @param value A list containing the base model as well as its parameters. It
#' must contains `model`, `par_base`, `cov_base`, `lag_base`, and `horizon`.
#' @export
`base<-` <- function(x, value) {

    base_res <- list(
        par_base = value$par_base,
        cov_base = value$cov_base,
        lag_base = value$lag_base
    )

    attr(x, "base") <- fit_base$model
    attr(x, "base_res") <- base_res
    attr(x, "horizon") <- value$horizon
    return(x)
}
