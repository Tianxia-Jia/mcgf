#' Compute the objective for wls method
#'
#' @param par Parameters of `cor_fn`.
#' @param cor_fn Correlation function.
#' @param cor_emp Empirical correlations.
#' @param par_fixed Fixed parameters of `cor_fn`.
#'
#' @keywords internal
#' @return The objective of weighted least squares.
obj_wls <- function(par, cor_fn, cor_emp, par_fixed) {

    fitted <- do.call(cor_fn, c(par, par_fixed))
    summand <- ((cor_emp - fitted) / (1 - fitted)) ^ 2
    summand[is.infinite(summand)] <- NA
    wls <- sum(summand, na.rm = T)
    return(wls)
}

#' Optimization for wls method
#'
#' @param par_init Initial values for parameters to be optimized.
#' @param method Parameter estimation method. "wls" or "mle".
#' @param optim_fn Optimization function.
#' @param cor_fn Correlation function.
#' @param par_fixed Fixed parameters of `cor_fn`.
#' @param lower Lower bound.
#' @param upper Upper bound.
#' @param ... Additional arguments passed to `optim_fn`, `obj_wls` or `obj_mle`.
#'
#' @keywords internal
#' @return A list outputted from optimization functions of `optim_fn`.
estimate <- function(par_init, method, optim_fn, cor_fn, par_fixed,
                     lower, upper, ...) {

    obj_fn <- switch(method,
                     wls = obj_wls,
                     mle = obj_mle)

    args <- list(par_init,
                 obj_fn,
                 cor_fn = cor_fn,
                 par_fixed = par_fixed,
                 lower = lower,
                 upper = upper,
                 ...)

    if (optim_fn == "optim") args <- c(args, method = "L-BFGS-B")

    return(do.call(optim_fn, args))
}
