#' Fit correlation base models
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of estimated parameters.
#' @export
#'
#' @details
#' Refer to [`fit_base.mcgf()`] and [`fit_base.mcgf_rs()`] for more details.
fit_base <- function(x, ...) {
    UseMethod("fit_base")
}

#' Parameter estimation for symmetric correlation functions for an `mcgf`
#' object.
#'
#' @param x An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
#' `sds`.
#' @param lag Integer time lag.
#' @param horizon Integer forecast horizon.
#' @param model Base model, one of `spatial`, `temporal`, `sep`, `fs`, `none`.
#' Only `sep` and `fs` are supported when `method = mle`. If `none`, NULLs are
#' returned.
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
#' @param dists_base List of distance matrices. If NULL, `dists(x)` is used.
#' Must be a matrix or an array of distance matrices.
#' @param scale_time Scale of time unit, default is 1. `lag` is divided by
#' `scale_time` for parameter estimation.
#' @param ... Additional arguments passed to `optim_fn`.
#'
#' @return A list containing outputs from optimization functions of `optim_fn`.
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
#' Note that both `wls` and `mle` are heuristic approaches when `x` contains
#' observations from a subset of the discrete spatial domain, though estimation
#' results are close to that using the full spatial domain for large sample
#' sizes.
#'
#' @examples
#' data(sim1)
#' sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#' sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
#' sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = 5)
#'
#' # Fit a pure spatial model
#' fit_spatial <- fit_base(
#'     sim1_mcgf,
#'     model = "spatial",
#'     lag = 5,
#'     par_init = c(c = 0.001, gamma = 0.5),
#'     par_fixed = c(nugget = 0)
#' )
#' fit_spatial$fit
#'
#' # Fit a pure temporal model
#' fit_temporal <- fit_base(
#'     sim1_mcgf,
#'     model = "temporal",
#'     lag = 5,
#'     par_init = c(a = 0.3, alpha = 0.5)
#' )
#' fit_temporal$fit
#'
#' # Fit a separable model
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
#' fit_sep$fit
#' @family functions related to fitting an `mcgf` object
fit_base.mcgf <- function(x,
                          lag,
                          horizon = 1,
                          model = c("spatial", "temporal", "sep", "fs", "none"),
                          method = c("wls", "mle"),
                          optim_fn = c("nlminb", "optim", "other"),
                          par_fixed = NULL,
                          par_init,
                          lower = NULL,
                          upper = NULL,
                          other_optim_fn = NULL,
                          dists_base = NULL,
                          scale_time = 1,
                          ...) {
    scale_time <- as.integer(scale_time)
    model <- match.arg(model)

    if (model == "none") {
        return(list(
            model = model,
            method = NULL,
            optim_fn = NULL,
            lag = NULL,
            fit = NULL,
            par_names = NULL,
            par_fixed = NULL,
            dists_base = NULL,
            scale_time = 1,
            dots = NULL
        ))
    }

    method <- match.arg(method)
    dots <- list(...)

    if (!is_numeric_scalar(lag)) {
        stop("`lag` must be a positive number.", call. = FALSE)
    } else if (lag < 0) {
        stop("`lag` must be a positive number.", call. = FALSE)
    }

    if (!is_numeric_scalar(horizon)) {
        stop("`horizon` must be a positive number.", call. = FALSE)
    } else if (horizon < 0) {
        stop("`horizon` must be a positive number.", call. = FALSE)
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

    lag_max <- lag + horizon - 1

    if (lag_max + 1 > length(acfs(x))) {
        stop("`lag` + `horizon` must be no greater than ", length(acfs(x)),
            ", or recompute `acfs` and `ccfs` with greater `lag_max`.",
            call. = FALSE
        )
    }

    if (method == "mle" && model %in% c("spatial", "temporal")) {
        stop("mle is available for `sep` and `fs` models only.", call. = FALSE)
    }

    par_model <- eval(as.name(paste0("par_", model)))
    lower_model <- eval(as.name(paste0("lower_", model)))
    upper_model <- eval(as.name(paste0("upper_", model)))

    if (!is.null(par_fixed)) {
        par_fixed_nm <- names(par_fixed)
        if (any(!par_fixed_nm %in% par_model)) {
            stop("unknow parameters in `par_fixed`.", call. = FALSE)
        }

        ind_not_fixed <- which(!par_model %in% par_fixed_nm)
        par_model <- par_model[ind_not_fixed]
        if (!is.null(lower)) {
            if (length(lower) != length(par_model)) {
                stop("`lower` must be of length ", length(par_model), ".",
                    call. = FALSE
                )
            }
            lower_model <- lower
        } else {
            lower_model <- lower_model[ind_not_fixed]
        }

        if (!is.null(upper)) {
            if (length(upper) != length(par_model)) {
                stop("`upper` must be of length ", length(par_model), ".",
                    call. = FALSE
                )
            }
            upper_model <- upper
        } else {
            upper_model <- upper_model[ind_not_fixed]
        }
    } else {
        par_fixed <- NULL

        if (!is.null(lower)) {
            if (length(lower) != length(par_model)) {
                stop("`lower` must be of length ", length(par_model), ".",
                    call. = FALSE
                )
            }
            lower_model <- lower
        }

        if (!is.null(upper)) {
            if (length(upper) != length(par_model)) {
                stop("`upper` must be of length ", length(par_model), ".",
                    call. = FALSE
                )
            }
            upper_model <- upper
        }
    }

    if (missing(par_init)) {
        stop("must provide `par_init`.", call. = FALSE)
    }

    par_init_nm <- names(par_init)
    if (any(!par_init_nm %in% par_model)) {
        stop("unknow parameters in `par_init`.", call. = FALSE)
    }

    if (any(!par_model %in% par_init_nm)) {
        par_missing <- par_model[which(!par_model %in% par_init_nm)]
        stop("initial value(s) for ",
            paste0("`", par_missing, "`", collapse = ", "),
            " not found.",
            call. = FALSE
        )
    }
    par_init <- par_init[order(match(names(par_init), par_model))]

    optim_fn <- match.arg(optim_fn)
    if (optim_fn == "other") {
        if (is.null(other_optim_fn)) {
            stop("specify a optimization function.", call. = FALSE)
        }
        optim_fn <- other_optim_fn
    }

    if (is.null(dists_base)) {
        dists_h <- dists(x)$h
    } else {
        if (model != "temporal") {
            check_dist(dists_base)

            if (is.matrix(dists_base)) {
                dists_h <- dists_base
            } else {
                if (model != "spatial") {
                    if (dim(dists_base)[3] < lag_max + 1) {
                        stop("third dim in `dists_base` must be greater or ",
                            "equal to ",
                            lag_max + 1,
                            ".",
                            call. = FALSE
                        )
                    }
                    dists_h <- dists_base[, , 1:(lag_max + 1)]
                } else {
                    dists_h <- dists_base[, , 1]
                }
            }
        }
    }

    if (method == "wls") {
        model_args <- switch(model,
            spatial = {
                cor_fn <- ".cor_exp"
                cor_emp <- ccfs(x)[, , 1]
                par_fixed_other <- list(x = dists_h)
                list(
                    cor_fn = cor_fn,
                    cor_emp = cor_emp,
                    par_fixed_other = par_fixed_other
                )
            },
            temporal = {
                cor_fn <- ".cor_cauchy"
                cor_emp <- acfs(x)[1:(lag_max + 1)]
                par_fixed_other <- list(x = 0:lag_max / scale_time)
                list(
                    cor_fn = cor_fn,
                    cor_emp = cor_emp,
                    par_fixed_other = par_fixed_other
                )
            },
            sep = {
                cor_fn <- "..cor_sep"
                cor_emp <- ccfs(x)[, , 1:(lag_max + 1)]
                h_u_ar <-
                    to_ar(h = dists_h, lag_max = lag_max)
                par_fixed_other <-
                    list(
                        h = h_u_ar$h_ar,
                        u = h_u_ar$u_ar / scale_time
                    )
                list(
                    cor_fn = cor_fn,
                    cor_emp = cor_emp,
                    par_fixed_other = par_fixed_other
                )
            },
            fs = {
                cor_fn <- ".cor_fs"
                cor_emp <- ccfs(x)[, , 1:(lag_max + 1)]
                h_u_ar <-
                    to_ar(h = dists_h, lag_max = lag_max)
                par_fixed_other <-
                    list(
                        h = h_u_ar$h_ar,
                        u = h_u_ar$u_ar / scale_time
                    )
                list(
                    cor_fn = cor_fn,
                    cor_emp = cor_emp,
                    par_fixed_other = par_fixed_other
                )
            }
        )

        res_base <- estimate(
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
    } else {
        model_args <- switch(model,
            sep = {
                cor_fn <- "..cor_sep"
                h_u_ar <-
                    to_ar(h = dists_h, lag_max = lag_max)
                par_fixed_other <-
                    list(
                        h = h_u_ar$h_ar,
                        u = h_u_ar$u_ar / scale_time
                    )
                list(
                    cor_fn = cor_fn,
                    par_fixed_other = par_fixed_other
                )
            },
            fs = {
                cor_fn <- ".cor_fs"
                h_u_ar <-
                    to_ar(h = dists_h, lag_max = lag_max)
                par_fixed_other <-
                    list(
                        h = h_u_ar$h_ar,
                        u = h_u_ar$u_ar / scale_time
                    )
                list(
                    cor_fn = cor_fn,
                    par_fixed_other = par_fixed_other
                )
            }
        )

        res_base <- estimate(
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
    }
    return(list(
        model = model,
        method = method,
        optim_fn = optim_fn,
        lag = lag,
        horizon = horizon,
        fit = res_base,
        par_names = names(par_init),
        par_fixed = par_fixed,
        dists_base = dists_base,
        scale_time = scale_time,
        dots = dots
    ))
}

#' Parameter estimation for symmetric correlation functions for an `mcgf_rs`
#' object.
#'
#' @param x An `mcgf_rs` object containing attributes `dists`, `acfs`, `ccfs`,
#' and `sds`.
#' @param lag_ls List of integer time lags.
#' @param horizon Integer forecast horizon.
#' @param model_ls List of base models, each element must be one of `spatial`,
#' `temporal`, `sep`, `fs`. Only `sep` and `fs` are supported when `mle` is used
#' in `model_ls`.
#' @param method_ls List of parameter estimation methods, weighted least square
#' (`wls`) or maximum likelihood estimation (`mle`).
#' @param optim_fn_ls List of optimization functions, each element must be one
#' of `nlminb`, `optim`, `other`. When use `other`, supply `other_optim_fn_ls`.
#' @param par_fixed_ls List of fixed parameters.
#' @param par_init_ls List of initial values for parameters to be optimized.
#' @param lower_ls Optional; list of lower bounds of parameters.
#' @param upper_ls Optional: list of upper bounds of parameters.
#' @param other_optim_fn_ls Optional, list of other optimization functions. The
#' first two arguments must be initial values for the parameters and a function
#' to be minimized respectively (same as that of `optim` and `nlminb`).
#' @param dists_base_ls List of lists of distance matrices. If NULL, `dists(x)`
#' is used. Each element must be a matrix or an array of distance matrices.
#' @param rs Logical; if TRUE `x` is treated as a regime-switching, FALSE if the
#' parameters need to be estimated in a non-regime-switching setting.
#' @param scale_time Scale of time unit, default is 1. `lag` is divided by
#' `scale_time` for parameter estimation.
#' @param ... Additional arguments passed to all `optim_fn_ls`.
#'
#' @return A list containing outputs from optimization functions of `optim_fn`
#' for each regime.
#' @export
#'
#' @details
#' This functions is the regime-switching variant of [`fit_base.mcgf()`].
#' Arguments are in lists. The length of arguments that end in `_ls` must be 1
#' or the same as the number of regimes in `x`. If the length of an argument is
#' 1, then it is set the same for all regimes. Refer to [`fit_base.mcgf()`] for
#' more details of the arguments.
#'
#' Note that both `wls` and `mle` are heuristic approaches when `x` contains
#' observations from a subset of the discrete spatial domain, though estimation
#' results are close to that using the full spatial domain for large sample
#' sizes.
#'
#' @examples
#' data(sim2)
#' sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#' sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
#' sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)
#'
#' # Fit a regime-switching pure spatial model
#' fit_spatial <-
#'     fit_base(
#'         sim2_mcgf,
#'         lag_ls = 5,
#'         model_ls = "spatial",
#'         par_init_ls = list(c(c = 0.000001, gamma = 0.5)),
#'         par_fixed_ls = list(c(nugget = 0))
#'     )
#' lapply(fit_spatial[1:2], function(x) x$fit)
#'
#' # Fit a regime-switching pure temporal model
#' fit_temporal <-
#'     fit_base(
#'         sim2_mcgf,
#'         lag_ls = 5,
#'         model_ls = "temporal",
#'         par_init_ls = list(
#'             list(a = 0.8, alpha = 0.8),
#'             list(a = 0.5, alpha = 0.5)
#'         )
#'     )
#' lapply(fit_temporal[1:2], function(x) x$fit)
#'
#' # Fit a regime-switching separable model
#' fit_sep <- fit_base(
#'     sim2_mcgf,
#'     lag_ls = 5,
#'     model_ls = "sep",
#'     par_init_ls = list(list(
#'         c = 0.000001,
#'         gamma = 0.5,
#'         a = 0.5,
#'         alpha = 0.5
#'     )),
#'     par_fixed_ls = list(c(nugget = 0))
#' )
#' lapply(fit_sep[1:2], function(x) x$fit)
#' @family functions related to fitting an `mcgf_rs` object
fit_base.mcgf_rs <- function(x,
                             lag_ls,
                             horizon = 1,
                             model_ls,
                             method_ls = "wls",
                             optim_fn_ls = "nlminb",
                             par_fixed_ls = list(NULL),
                             par_init_ls,
                             lower_ls = list(NULL),
                             upper_ls = list(NULL),
                             other_optim_fn_ls = list(NULL),
                             dists_base_ls = list(NULL),
                             scale_time = 1,
                             rs = TRUE,
                             ...) {
    scale_time <- as.integer(scale_time)

    args_ls <- c(
        "lag", "model", "method", "optim_fn", "par_fixed", "par_init",
        "lower", "upper", "other_optim_fn", "dists_base"
    )
    args_i <- paste0("i_", args_ls)
    args_rs <- paste0(args_ls, "_ls")

    if (rs) {
        lvs <- levels((attr(x, "label", exact = TRUE)))
        n_regime <- length(lvs)
        res_base_ls <- vector("list", n_regime)

        for (i in 1:length(args_rs)) {
            length_args_i <- length(eval(as.name(args_rs[i])))

            if (length_args_i == 1) {
                assign(args_i[i], rep(1L, n_regime))
            } else if (length_args_i == n_regime) {
                assign(args_i[i], 1:n_regime)
            } else {
                stop("length of `", args_rs[i], "` must be 1 or ", n_regime,
                    ".",
                    call. = FALSE
                )
            }
        }

        fit_base_fixed <- list(horizon = horizon, scale_time = scale_time, ...)

        for (n in 1:n_regime) {
            ind_n <- lapply(mget(args_i), function(x) x[[n]])
            args_rs_n <- mget(args_rs)
            names(args_rs_n) <- args_ls
            args_n <- Map(function(x, ind) x[[ind]], args_rs_n, ind_n)

            x_n <- x
            acfs(x_n) <- acfs(x)$acfs_rs[[n]]
            ccfs(x_n) <- ccfs(x)$ccfs_rs[[n]]
            sds(x_n) <- sds(x)$sds_rs[[n]]
            attr(x_n, "mle_label") <- lvs[n]
            fit_base_fixed_n <- c(
                fit_base_fixed,
                list(x = x_n)
            )
            res_base_ls[[n]] <- do.call(
                fit_base.mcgf,
                c(args_n, fit_base_fixed_n)
            )
        }
        names(res_base_ls) <- paste0("Regime ", lvs)
        res_base_ls <- c(res_base_ls, rs = rs)
        return(res_base_ls)
    } else {
        for (i in 1:length(args_rs)) {
            value_args_i <- eval(as.name(args_rs[i]))[[1]]
            assign(args_ls[i], value_args_i)
        }

        args_no_rs <- mget(args_ls)
        names(args_no_rs) <- args_ls

        x_no_rs <- x
        acfs(x_no_rs) <- acfs(x)$acfs
        ccfs(x_no_rs) <- ccfs(x)$ccfs
        sds(x_no_rs) <- sds(x)$sds
        fit_base_fixed <- c(
            horizon = horizon,
            list(x = x_no_rs, scale_time = scale_time),
            ...
        )
        res_base_ls <- do.call(
            fit_base.mcgf,
            c(args_no_rs, fit_base_fixed)
        )
        res_base_ls <- c(list(res_base_ls), rs = rs)
        return(res_base_ls)
    }
}
