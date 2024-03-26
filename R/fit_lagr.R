#' Fit correlation Lagrangian models
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return A vector of estimated parameters.
#' @export
#'
#' @details
#' Refer to [`fit_lagr.mcgf()`] and [`fit_lagr.mcgf_rs()`] for more details.
fit_lagr <- function(x, ...) {
    UseMethod("fit_lagr")
}

#' Parameter estimation for Lagrangian correlation functions for an `mcgf`
#' object.
#'
#' @param x An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
#' `sds`. `x` must have been passed to `add_base()` or `base<-`
#' @param model Base model, one of `lagr_tri`, `lagr_askey`, `lagr_exp`, or
#' `none`.
#' If `none`, NULLs are returned.
#' @param method Parameter estimation methods, weighted least square (`wls`) or
#' maximum likelihood estimation (`mle`).
#' @param optim_fn Optimization functions, one of `nlminb`, `optim`, `other`.
#' When `optim_fn = other`, supply `other_optim_fn`.
#' @param par_fixed Fixed parameters.
#' @param par_init Initial values for parameters to be optimized.
#' @param lower Optional; lower bounds of parameters lambda, v1, v2, and k.
#' @param upper Optional: upper bounds of parameters lambda, v1, v2, and k.
#' @param other_optim_fn Optional, other optimization functions. The first two
#' arguments must be initial values for the parameters and a function to be
#' minimized respectively (same as that of `optim` and `nlminb`).
#' @param dists_base Logical; if TRUE `dists_base` from the base model is used
#' as the distance.
#' @param dists_lagr List of distance matrices/arrays. Used when `dists_base` is
#' FALSE. If NULL, `dists(x)` is used.
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
#' Note that both `wls` and `mle` are heuristic approaches when `x` contains
#' observations from a subset of the discrete spatial domain, though estimation
#' results are close to that using the full spatial domain for large sample
#' sizes.
#'
#' Since parameters for the base model and the Lagrangian model are estimated
#' sequentially, more accurate estimation may be obtained if the full model is
#' fitted all at once.
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
#' fit_lagr$fit
#' @family functions on fitting an mcgf
fit_lagr.mcgf <- function(x,
                          model = c(
                              "lagr_tri", "lagr_askey", "lagr_exp",
                              "none"
                          ),
                          method = c("wls", "mle"),
                          optim_fn = c("nlminb", "optim", "other"),
                          par_fixed = NULL,
                          par_init,
                          lower = NULL,
                          upper = NULL,
                          other_optim_fn = NULL,
                          dists_base = FALSE,
                          dists_lagr = NULL,
                          ...) {
    model <- match.arg(model)

    if (model == "none") {
        return(list(
            model = model,
            method = NULL,
            optim_fn = NULL,
            fit = NULL,
            par_names = NULL,
            par_fixed = NULL,
            dists_base = NULL,
            dists_lagr = NULL,
            dots = NULL
        ))
    }

    if (!is.null(lower)) lower <- unlist(lower)
    if (!is.null(upper)) upper <- unlist(upper)

    method <- match.arg(method)
    dots <- list(...)

    par_model <- c("lambda", "v1", "v2", "k")
    lower_model <- c(0, -99999, -99999, 0)
    upper_model <- c(1, 99999, 99999, 99999)

    lag <- attr(x, "lag", exact = TRUE)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_max <- lag + horizon - 1
    scale_time <- attr(x, "scale_time", exact = TRUE)

    if (!is.null(par_fixed)) {
        par_fixed_nm <- names(par_fixed)

        if (is.null(par_fixed_nm) || any(!par_fixed_nm %in% par_model)) {
            stop("unknow parameters in `par_fixed`.")
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

    if (is.null(par_init)) {
        stop("must provide `par_init`.", call. = FALSE)
    }

    par_init_nm <- names(par_init)
    if (is.null(par_init_nm) || any(!par_init_nm %in% par_model)) {
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

    if (dists_base) {
        if (!is.null(dists_lagr)) {
            warning("`dists_base` is used. Set it to FALSE if `dists_lagr` ",
                "needs to be used.",
                call. = FALSE
            )
        }
        dists_base_ls <- attr(x, "base_res", exact = TRUE)$dists_base

        if (is.null(dists_base_ls)) {
            stop("`dists_base` is NULL.", call. = FALSE)
        }

        if (any(names(dists_base_ls) != c("h", "h1", "h2"))) {
            stop("`dists_base` must be a list of matrices/arrays with ",
                "names `h`, `h1`, 'h2'.",
                call. = FALSE
            )
        }

        lagr_h <- dists_base_ls$h
        lagr_h1 <- dists_base_ls$h1
        lagr_h2 <- dists_base_ls$h2
    } else {
        if (is.null(dists_lagr)) {
            lagr_h <- dists(x)$h
            lagr_h1 <- dists(x)$h1
            lagr_h2 <- dists(x)$h2
        } else {
            if (any(names(dists_lagr) != c("h", "h1", "h2"))) {
                stop("`dists_lagr` must be a list of matrices/arrays with ",
                    "names `h`, `h1`, 'h2'.",
                    call. = FALSE
                )
            } else {
                check_dist(dists_lagr$h, "h", check_sym = FALSE)
                check_dist_sign(dists_lagr$h1, "h1", check_sym = FALSE)
                check_dist_sign(dists_lagr$h2, "h2", check_sym = FALSE)

                lagr_h <- dists_lagr$h
                lagr_h1 <- dists_lagr$h1
                lagr_h2 <- dists_lagr$h2
            }
        }
    }

    if (any(dim(lagr_h) != dim(lagr_h1))) {
        stop("unmatching dimensions for `h` and `h1 `in `dists`.",
            call. = FALSE
        )
    }
    if (any(dim(lagr_h) != dim(lagr_h2))) {
        stop("unmatching dimensions for `h` and `h2 `in `dists`.",
            call. = FALSE
        )
    }

    if (!is.matrix(lagr_h)) {
        if (dim(lagr_h)[3] < lag_max + 1) {
            stop("third dims from `dists_lagr` must be greater or ",
                "equal to ", lag_max + 1, ".",
                call. = FALSE
            )
        }

        lagr_h <- lagr_h[, , 1:(lag_max + 1)]
        lagr_h1 <- lagr_h1[, , 1:(lag_max + 1)]
        lagr_h2 <- lagr_h2[, , 1:(lag_max + 1)]
    }

    cor_base <- attr(x, "base_res", exact = TRUE)$cor_base
    u_ar <- to_ar(h = lagr_h, lag_max = lag_max)$u_ar
    h1_ar <- to_ar(h = lagr_h1, lag_max = lag_max, u = FALSE)
    h2_ar <- to_ar(h = lagr_h2, lag_max = lag_max, u = FALSE)

    par_fixed_other <- list(
        cor_base = cor_base,
        lagrangian = model,
        h1 = h1_ar,
        h2 = h2_ar,
        u = u_ar / scale_time
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
            lag = lag,
            ...
        )
    }

    return(list(
        model = model,
        method = method,
        optim_fn = optim_fn,
        fit = res_lagr,
        par_names = names(par_init),
        par_fixed = par_fixed,
        dists_base = dists_base,
        dists_lagr = dists_lagr,
        dots = dots
    ))
}

#' Parameter estimation for Lagrangian correlation functions for an `mcgf_rs`
#' object.
#'
#' @param x An `mcgf_rs` object containing attributes `dists`, `acfs`, `ccfs`,
#' and `sds`. `x` must have been passed to `add_base()` or `base<-`
#' @param model_ls List of base models, each element must be one of `lagr_tri`,
#' `lagr_askey`, `lagr_exp`, or `none`. If `none`, NULLs are returned.
#' @param method_ls List of parameter estimation methods, weighted least square
#' (`wls`) or maximum likelihood estimation (`mle`).
#' @param optim_fn_ls List of optimization functions, each element must be one
#' of `nlminb`, `optim`, `other`. When use `other`, supply `other_optim_fn_ls`
#' @param par_fixed_ls List of fixed parameters.
#' @param par_init_ls List of initial values for parameters to be optimized.
#' @param lower_ls Optional; list of lower bounds of parameters.
#' @param upper_ls Optional: list of upper bounds of parameters.
#' @param other_optim_fn_ls Optional, list of other optimization functions. The
#' first two arguments must be initial values for the parameters and a function
#' to be minimized respectively (same as that of `optim` and `nlminb`).
#' @param dists_base_ls List of lists of distance matrices. If NULL, `dists(x)`
#' is used. Each element must be a matrix or an array of distance matrices.
#' @param dists_lagr_ls List of distance matrices/arrays. Used when
#' `dists_base` is FALSE. If NULL, `dists(x)` is used.
#' @param rs Logical; if TRUE `x` is treated as a regime-switching, FALSE if the
#' parameters need to be estimated in a non-regime-switching setting.
#' @param ... Additional arguments passed to `optim_fn`.
#'
#' @return A list containing outputs from optimization functions of `optim_fn`.
#' @export
#'
#' @details
#' This functions is the regime-switching variant of [`fit_lagr.mcgf()`].
#' Arguments are in lists. The length of arguments that end in `_ls` must be 1
#' or the same as the number of regimes in `x`. If the length of an argument is
#' 1, then it is set the same for all regimes. Refer to [`fit_lagr.mcgf()`] for
#' more details of the arguments.
#'
#' Note that both `wls` and `mle` are heuristic approaches when `x` contains
#' observations from a subset of the discrete spatial domain, though estimation
#' results are close to that using the full spatial domain for large sample
#' sizes.
#'
#' Since parameters for the base model and the Lagrangian model are estimated
#' sequentially, more accurate estimation may be obtained if the full model is
#' fitted all at once.
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
#' lapply(fit_lagr_rs[1:2], function(x) x$fit)
#' @family functions on fitting an mcgf_rs
fit_lagr.mcgf_rs <- function(x,
                             model_ls,
                             method_ls = "wls",
                             optim_fn_ls = "nlminb",
                             par_fixed_ls = list(NULL),
                             par_init_ls,
                             lower_ls = list(NULL),
                             upper_ls = list(NULL),
                             other_optim_fn_ls = list(NULL),
                             dists_base_ls,
                             dists_lagr_ls = list(NULL),
                             rs = TRUE,
                             ...) {
    if (missing(dists_base_ls)) {
        dists_base_ls <- list(FALSE)
    }

    args_ls <- c(
        "model", "method", "optim_fn", "par_fixed", "par_init",
        "lower", "upper", "other_optim_fn", "dists_base", "dists_lagr"
    )
    args_i <- paste0("i_", args_ls)
    args_rs <- paste0(args_ls, "_ls")

    lag_ls <- attr(x, "lag", exact = TRUE)

    base_rs <- attr(x, "base_rs", exact = TRUE)
    if (length(base_rs) != 1) base_rs <- any(base_rs)

    if (rs) {
        lvs <- levels((attr(x, "label", exact = TRUE)))
        n_regime <- length(lvs)
        res_lagr_ls <- vector("list", n_regime)

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

        if (length(lag_ls) == 1) {
            lag_ls <- rep(lag_ls, n_regime)
        }

        for (n in 1:n_regime) {
            ind_n <- lapply(mget(args_i), function(x) x[[n]])
            args_rs_n <- mget(args_rs)
            names(args_rs_n) <- args_ls
            args_n <- Map(function(x, ind) x[[ind]], args_rs_n, ind_n)

            x_n <- x
            acfs(x_n) <- acfs(x)$acfs_rs[[n]]
            ccfs(x_n) <- ccfs(x)$ccfs_rs[[n]]
            sds(x_n) <- sds(x)$sds_rs[[n]]
            attr(x_n, "lag") <- lag_ls[[n]]
            attr(x_n, "mle_label") <- lvs[n]
            if (base_rs) {
                attr(x_n, "base_res") <-
                    attr(x_n, "base_res", exact = TRUE)[[n]]
            } else {
                attr(x_n, "base_res") <- attr(x_n, "base_res", exact = TRUE)
            }
            res_lagr_ls[[n]] <- do.call(
                fit_lagr.mcgf,
                c(args_n, list(x = x_n), ...)
            )
        }

        names(res_lagr_ls) <- paste0("Regime ", lvs)
        res_lagr_ls <- c(res_lagr_ls, rs = rs)
        return(res_lagr_ls)
    } else {
        if (base_rs) {
            stop("the base model cannot be regime-switching if the Lagragian",
                " model is not regime-switching.",
                call. = FALSE
            )
        }

        if (length(unique(lag_ls)) != 1) {
            stop("`lag` must be the same for all regimes.", call. = FALSE)
        }

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
        attr(x_no_rs, "lag") <- lag_ls[[1]]
        if (base_rs) {
            attr(x_no_rs, "base_res") <-
                attr(x_no_rs, "base_res", exact = TRUE)[[1]]
        } else {
            attr(x_no_rs, "base_res") <- attr(x_no_rs, "base_res", exact = TRUE)
        }
        res_lagr_ls <- do.call(
            fit_lagr.mcgf,
            c(args_no_rs, list(x = x_no_rs), ...)
        )
        res_lagr_ls <- c(list(res_lagr_ls), rs = rs)
        return(res_lagr_ls)
    }
}
