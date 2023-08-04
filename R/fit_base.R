#' Fit correlation models
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of estimated parameters
#' @export
#' @family {functions related to model fitting}
fit_base <- function(x, ...) {
    UseMethod("fit_base")
}

par_spatial <- c("c", "gamma", "nugget")
lower_spatial <- c(0, 0, 0)
upper_spatial <- c(100, 0.5, 1)

par_temporal <- c("a", "alpha")
lower_temporal <- c(0, 0)
upper_temporal <- c(100, 1)

par_sep <- c(par_spatial, par_temporal)
lower_sep <- c(lower_spatial, lower_temporal)
upper_sep <- c(upper_spatial, upper_temporal)

par_fs <- c(par_sep, "beta")
lower_fs <- c(lower_sep, 0)
upper_fs <- c(upper_sep, 1)

#' Parameter estimation for symmetric correlation functions.
#'
#' @param x An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
#' `sds`.
#' @param lag Integer time lag.
#' @param horizon Integer forecast horizon.
#' @param model Base model, one of `spatial`, `temporal`, `sep`, `fs`. Only
#' `sep` and `fs` are supported when `method = mle`
#' @param method Parameter estimation methods, weighted least square (`wls`) or
#' maximum likelihood estimation (`mle`).
#' @param optim_fn Optimization functions, one of `nlminb`, `optim`, `other`.
#' When `optim_fn = other`, supply `other_optim_fn`.
#' @param par_fixed Fixed parameters.
#' @param par_init Initial values for parameters to be optimized.
#' @param lower Optional; lower bound of parameters.
#' @param upper Optional: upper bound of parameters.
#' @param other_optim_fn Optional, other optimization functions. The first two
#' arguments must be initial values for the parameters and a function to be
#' minimized respectively (same as that of `optim` and `nlminb`).
#' @param ... Additional arguments passed to `optim_fn`.
#'
#' @return A list outputted from optimization functions of `optim_fn`.
#' @export
#'
#' @details
#' This function fits the separable and fully symmetric models using weighted
#' least squares and maximum likelihood estimation. Optimization functions such
#' as `nlminb` and `optim` are supported. Other functions are also supported by
#' setting `optim_fn = "other"` and supplying `other_optim_fn`. `lower` and
#' `upper` are lower and upper bounds of parameters in `par_init` and default
#' bounds are used if they are not specified.
#'
#' @family {functions related to model fitting}
fit_base.mcgf <- function(x,
                          lag,
                          horizon = 1,
                          model = c("spatial", "temporal", "sep", "fs"),
                          method = c("wls", "mle"),
                          optim_fn = c("nlminb", "optim", "other"),
                          par_fixed,
                          par_init,
                          lower,
                          upper,
                          other_optim_fn,
                          ...) {

    if(!is_numeric_scalar(lag)) {
        stop("`lag` must be a positive number.")
    } else if (lag < 0) {
        stop("`lag` must be a positive number.")
    }

    if(!is_numeric_scalar(horizon)) {
        stop("`horizon` must be a positive number.")
    } else if (horizon < 0) {
        stop("`horizon` must be a positive number.")
    }

    lag_max <- lag + horizon - 1

    if (lag_max + 1 > length(acfs(x)))
        stop("`lag` + `horizon` must be no greater than ", length(acfs(x)),
             ", or recompute `acfs` and `ccfs` with greater `lag_max`.")

    model <- match.arg(model)
    method <- match.arg(method)

    if (method == "mle" && model %in% c("spatial", "temporal"))
        stop("mle is available for `sep` and `fs` models only.")

    par_model <- eval(as.name(paste0("par_", model)))
    lower_model <- eval(as.name(paste0("lower_", model)))
    upper_model <- eval(as.name(paste0("upper_", model)))

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

    if (method == "wls") {

        model_args <- switch(
            model,
            spatial = {
                cor_fn <- ".cor_exp"
                cor_emp <- ccfs(x)[, , 1]
                par_fixed_other <- list(x = dists(x)$h)
                list(cor_fn = cor_fn,
                     cor_emp = cor_emp,
                     par_fixed_other = par_fixed_other)

            },
            temporal = {
                cor_fn <- ".cor_cauchy"
                cor_emp <- acfs(x)[1:(lag_max + 1)]
                par_fixed_other <<- list(x = 0:lag_max, nu = 1, nugget = 0)
                list(cor_fn = cor_fn,
                     cor_emp = cor_emp,
                     par_fixed_other = par_fixed_other)
            },
            sep = {
                cor_fn <- "..cor_sep"
                cor_emp <- ccfs(x)[, , 1:(lag_max + 1)]
                h_u_ar <-
                    to_ar(h = dists(x)$h, lag_max = lag_max)
                par_fixed_other <-
                    list(h = h_u_ar$h_ar,
                         u = h_u_ar$u_ar)
                list(cor_fn = cor_fn,
                     cor_emp = cor_emp,
                     par_fixed_other = par_fixed_other)
            },
            fs = {
                cor_fn <- ".cor_fs"
                cor_emp <- ccfs(x)[, , 1:(lag_max + 1)]
                h_u_ar <-
                    to_ar(h = dists(x)$h, lag_max = lag_max)
                par_fixed_other <-
                    list(h = h_u_ar$h_ar,
                         u = h_u_ar$u_ar)
                list(cor_fn = cor_fn,
                     cor_emp = cor_emp,
                     par_fixed_other = par_fixed_other)
            }
        )

        res_wls <- estimate(
            par_init = par_init,
            method = method,
            optim_fn = optim_fn,
            cor_fn = model_args$cor_fn,
            cor_emp = model_args$cor_emp,
            par_fixed = c(par_fixed, model_args$par_fixed_other),
            lower = lower_model,
            upper = upper_model,
            ...
        )
        return(c(res_wls, par_names = list(names(par_init)),
                 list(par_fixed = par_fixed)))

    } else {

        model_args <- switch(
            model,
            sep = {
                cor_fn <- "..cor_sep"
                h_u_ar <-
                    to_ar(h = dists(x)$h, lag_max = lag_max)
                par_fixed_other <-
                    list(h = h_u_ar$h_ar,
                         u = h_u_ar$u_ar)
                list(cor_fn = cor_fn,
                     par_fixed_other = par_fixed_other)
            },
            fs = {
                cor_fn <- ".cor_fs"
                h_u_ar <-
                    to_ar(h = dists(x)$h, lag_max = lag_max)
                par_fixed_other <-
                    list(h = h_u_ar$h_ar,
                         u = h_u_ar$u_ar)
                list(cor_fn = cor_fn,
                     par_fixed_other = par_fixed_other)
            }
        )

        res_mle <- estimate(
            par_init = par_init,
            method = method,
            optim_fn = optim_fn,
            cor_fn = model_args$cor_fn,
            par_fixed = c(par_fixed, model_args$par_fixed_other),
            lower = lower_model,
            upper = upper_model,
            x = x,
            lag = lag,
            ...
        )
        return(c(res_mle, par_names = list(names(par_init)),
                 list(par_fixed = par_fixed)))
    }
}
