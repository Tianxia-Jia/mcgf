#' Fit correlation models
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of estimated parameters
#' @export
#' @family {functions related to the class}
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
#' Additional details...
#'
#' @examples
#' wind_sq <- sqrt(wind[, -1])
#' time <- wind[, 1]
#' wind_mcgf <- mcgf(data = wind_sq, locations = wind_loc, time = time)
#' wind_mcgf <- add_acfs(x = wind_mcgf, lag_max = 3)
#' wind_mcgf <- add_ccfs(x = wind_mcgf, lag_max = 3)
#' fit_base(wind_mcgf, lag = 3, model = "spatial", optim_fn = "nlminb",
#'          par_init = list(nugget = 0.01, c = 0.001), par_fixed = list(gamma = 0.5))
#' fit_base(wind_mcgf, lag = 3, model = "spatial", optim_fn = "optim",
#'          par_init = list(nugget = 0.01, c = 0.001), par_fixed = list(gamma = 0.5))
#' fit_base(wind_mcgf, lag = 3, model = "temporal", optim_fn = "optim",
#'          par_init = list(a = 1, alpha = 0.8))
#' fit_base(wind_mcgf, lag = 3, model = "temporal", optim_fn = "nlminb",
#'          par_init = list(a = 1, alpha = 0.8))
#' fit_base(wind_mcgf, lag = 3, model = "sep", optim_fn = "nlminb",
#'          par_init = list(a = 1, alpha = 1, nugget = 0.02, c = 0.001),
#'          par_fixed = list(gamma = 0.5))
#' fit_base(wind_mcgf, lag = 3, model = "fs", optim_fn = "nlminb",
#'          par_init = list(beta = 0),
#'          par_fixed = list(a = 0.972, alpha = 0.834, nugget = 0.0415,
#'          c = 0.00128, gamma = 0.5))
#' fit_base(wind_mcgf, lag = 3, model = "fs", optim_fn = "nlminb",
#'          par_init = list(beta = 0),
#'          par_fixed = list(a = 0.808173954  , alpha = 0.818269388 ,
#'          nugget = 0.048333797  , c = 0.001110266  , gamma = 0.5))
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

    model <- match.arg(model)
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

    method <- match.arg(method)
    optim_fn <- match.arg(optim_fn)
    if (optim_fn == "other") {
        if (missing(other_optim_fn))
            stop("specify a optimization function.")
        optim_fn = other_optim_fn
    }

    if (method == "wls") {
        res <- optim_wls(par_init = par_init,
                         optim_fn = optim_fn,
                         cor_fn = model_args$cor_fn,
                         cor_emp = model_args$cor_emp,
                         par_fixed = c(par_fixed, model_args$par_fixed_other),
                         lower = lower_model,
                         upper = upper_model,
                         ...)
        return(res)
    }
}

#' Compute the objective for wls method
#'
#' @param par Parameters of `cor_fn`.
#' @param cor_fn Correlation function.
#' @param cor_emp Empirical correlations.
#' @param par_fixed Fixed parameters of `cor_fn`.
#'
#' @keywords internal
#' @return The objective of weighted least squares.
cor_wls <- function(par, cor_fn, cor_emp, par_fixed) {

    fitted <- do.call(cor_fn, c(par, par_fixed))
    summand <- ((cor_emp - fitted) / (1 - fitted)) ^ 2
    summand[is.infinite(summand)] <- NA
    wls <- sum(summand, na.rm = T)
    return(wls)
}

#' Optimization for wls method
#'
#' @param par_init Initial values for parameters to be optimized.
#' @param optim_fn Optimization function.
#' @param cor_fn Correlation function.
#' @param cor_emp Empirical correlations.
#' @param par_fixed Fixed parameters of `cor_fn`.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param ... Additional arguments passed to `optim_fn`.
#'
#' @keywords internal
#' @return A list outputted from optimization functions of `optim_fn`.
optim_wls <- function(par_init, optim_fn, cor_fn, cor_emp, par_fixed, lower,
                      upper, ...) {

    args <- list(par_init,
                 cor_wls,
                 cor_fn = cor_fn,
                 cor_emp = cor_emp,
                 par_fixed = par_fixed,
                 lower = lower,
                 upper = upper,
                 ...)

    if (optim_fn == "optim") args <- c(args, method = "L-BFGS-B")

    return(do.call(optim_fn, args))
}
