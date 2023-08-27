#' Fit correlation Lagrangian models
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of estimated parameters
#' @export
#' @family {functions related to model fitting}
fit_lagr <- function(x, ...) {
    UseMethod("fit_lagr")
}

#' Parameter estimation for Lagrangian correlation functions.
#'
#' @param x An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
#' `sds`. `x` must have been passed to `add_base()` or `base<-`
#' @param model Base model, one of `lagr_tri`, `lagr_askey`.
#' @param method Parameter estimation methods, weighted least square (`wls`) or
#' maximum likelihood estimation (`mle`).
#' @param optim_fn Optimization functions, one of `nlminb`, `optim`, `other`.
#' When `optim_fn = other`, supply `other_optim_fn`.
#' @param par_fixed Fixed parameters.
#' @param par_init Initial values for parameters to be optimized.
#' @param lower Optional; lower bounds of parameters.
#' @param upper Optional: upper bounds of parameters.
#' @param other_optim_fn Optional, other optimization functions. The first two
#' arguments must be initial values for the parameters and a function to be
#' minimized respectively (same as that of `optim` and `nlminb`).
#' @param ... Additional arguments passed to `optim_fn`.
#'
#' @return A list containing outputs from optimization functions of `optim_fn`.
#' @export
#'
#' @details
#' This function fits the Lagrangian models using weighted least squares and
#' maximum likelihood estimation. The base model must be fitted first using
#' `add_base()` or `base<-`. Optimization functions such as `nlminb` and `optim`
#' are supported. Other functions are also supported by setting
#' `optim_fn = "other"` and supplying `other_optim_fn`. `lower` and `upper` are
#' lower and upper bounds of parameters in `par_init` and default bounds are
#' used if they are not specified.
#'
#' @family {functions related to model fitting}
fit_lagr.mcgf <- function(x,
                          model = c("lagr_tri", "lagr_askey"),
                          method = c("wls", "mle"),
                          optim_fn = c("nlminb", "optim", "other"),
                          par_fixed,
                          par_init,
                          lower,
                          upper,
                          other_optim_fn,
                          ...) {

    dots <- list(...)

    par_model <- c("lambda", "v1", "v2", "k")
    lower_model <- c(0, -9999, -9999, 0)
    upper_model <- c(1, 9999, 9999, 9999)

    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_max <- lag + horizon - 1

    model <- match.arg(model)
    method <- match.arg(method)

    if (!missing(par_fixed)) {
        par_fixed_nm <- names(par_fixed)
        if (is.null(par_fixed_nm) || any(!par_fixed_nm %in% par_model))
            stop("unknow parameters in `par_fixed`.")

        ind_not_fixed <- which(!par_model %in% par_fixed_nm)
        par_model <- par_model[ind_not_fixed]
        if (!missing(lower)) {
            if (length(lower) != length(par_model))
                stop("`lower` must be of length ", length(par_model), ".")
            lower_model <- lower
        }
        lower_model <- lower_model[ind_not_fixed]

        if (!missing(upper)) {
            if (length(upper) != length(par_model))
                stop("`upper` must be of length ", length(par_model), ".")
            upper_model <- upper
        }
        upper_model <- upper_model[ind_not_fixed]
    } else {
        par_fixed <- NULL
    }

    if (missing(par_init))
        stop("must provide `par_init`.")

    par_init_nm <- names(par_init)
    if (is.null(par_init_nm) || any(!par_init_nm %in% par_model))
        stop("unknow parameters in `par_init`.")

    if (any(!par_model %in% par_init_nm)) {
        par_missing <- par_model[which(!par_model %in% par_init_nm)]
        stop("initial value(s) for ",
             paste0('`', par_init_nm, "`", collapse = ", "),
             " not found.")
    }

    par_init <- par_init[order(match(names(par_init), par_model))]

    optim_fn <- match.arg(optim_fn)
    if (optim_fn == "other") {
        if (missing(other_optim_fn))
            stop("specify a optimization function.")
        optim_fn = other_optim_fn
    }

    cor_base <- attr(x, "base_res", exact = TRUE)$cor_base
    u_ar <- to_ar(h = dists(x)$h, lag_max = lag_max)$u_ar
    h1_ar <- to_ar(h = dists(x)$h1, lag_max = lag_max, u = FALSE)
    h2_ar <- to_ar(h = dists(x)$h2, lag_max = lag_max, u = FALSE)

    par_fixed_other <- list(
        cor_base = cor_base,
        lagrangian = model,
        h1 = h1_ar,
        h2 = h2_ar,
        u = u_ar
    )

    if (method == "wls") {

        res_lagr <- estimate(
            par_init = par_init,
            method = method,
            optim_fn = optim_fn,
            cor_fn = "..cor_stat",
            cor_emp = ccfs(x)[, , 1:(lag_max + 1)],
            par_fixed = c(par_fixed, par_fixed_other),
            lower = lower_model,
            upper = upper_model,
            ...
        )

    } else {

        res_lagr <- estimate(
            par_init = par_init,
            method = method,
            optim_fn = optim_fn,
            cor_fn = "..cor_stat",
            par_fixed = c(par_fixed, par_fixed_other),
            lower = lower_model,
            upper = upper_model,
            x = x,
            lag = lag
        )
    }

    return(list(
        model = model,
        method = method,
        optim_fn = optim_fn,
        fit = res_lagr,
        par_names = names(par_init),
        par_fixed = par_fixed,
        dots = dots))
}
