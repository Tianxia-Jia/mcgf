#' Calculate general stationary correlation.
#'
#' @param n_regime Integer, number of regimes.
#' @param base_ls List of base model, `sep` or `fs` for now. Or list of
#' correlation matrices/arrays.
#' @param lagrangian_ls List of Lagrangian model, `lagr_tri` for now.
#' @param par_base_ls List of parameters for the base model, used only when
#' `base_fixed = FALSE`.
#' @param par_lagr_ls List of parameters for the Lagrangian model.
#' @param lambda_ls List of weight of the Lagrangian term,
#' \eqn{\lambda\in[0, 1]}.
#' @param h_ls List of Euclidean distance matrix or array,
#' used only when `base_fixed = FALSE`.
#' @param h1_ls List of horizontal distance matrix or array, same dimension as
#' `h_ls`.
#' @param h2_ls List of vertical distance matrix or array, same dimension as
#' `h_ls`.
#' @param u_ls List of time lag, same dimension as `h_ls`.
#' @param base_fixed Logical; if TRUE, `base_ls` is the list of correlation.
#'
#' @return Correlations for the general stationary model. Same dimension of
#' `base_ls` if `base_fixed = TRUE`.
#' @export
#'
#' @details It gives a list of general stationary correlation for `n_regime`
#' regimes. See [cor_stat] for the model details. Model parameters are lists of
#' length 1 or `n_regime`. When length is 1, same values are used for all
#' regimes. If `base_fixed = TRUE`, the base is a list of correlation and
#' `par_base_ls` and `h_ls` are not used.
#'
#' @examples
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' par_t <- list(a = 1, alpha = 0.5)
#' par_base <- list(par_s = par_s, par_t = par_t)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
#' h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
#' h <- sqrt(h1 ^ 2 + h2 ^ 2)
#' u <- matrix(1, nrow = 2, ncol = 2)
#' cor_stat_rs(
#'     n_regime = 2,
#'     par_base_ls = list(par_base),
#'     par_lagr_ls = list(par_lagr),
#'     h_ls = list(h),
#'     h1_ls = list(h1),
#'     h2_ls = list(h2),
#'     u_ls = list(u, u + 1),
#'     lambda_ls = list(0.8),
#'     base_ls = list("sep"),
#'     lagrangian_ls = list("lagr_tri")
#' )
#'
#' h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
#' h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
#' h <- sqrt(h1 ^ 2 + h2 ^ 2)
#' u <- array(rep(c(1, 2, 3), each = 4), dim = c(2, 2, 3))
#' fit_base <- cor_fs(nugget = 0.5, c = 0.01, gamma = 0.5, a = 1, alpha = 0.5,
#'                    beta = 0.0, h = h, u = u)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' cor_stat_rs(
#'     n_regime = 2,
#'     par_lagr_ls = list(par_lagr),
#'     h1_ls = list(h1),
#'     h2_ls = list(h2),
#'     u_ls = list(u, u + 1),
#'     lambda_ls = list(0.8),
#'     base_ls = list(fit_base),
#'     lagrangian = list("lagr_tri"),
#'     base_fixed = TRUE
#' )
#'
#' @seealso [cor_stat]
cor_stat_rs <- function(n_regime,
                        base_ls,
                        lagrangian_ls,
                        par_base_ls,
                        par_lagr_ls,
                        lambda_ls,
                        h_ls,
                        h1_ls,
                        h2_ls,
                        u_ls,
                        base_fixed = FALSE) {
    args_stat <- names(formals(cor_stat))
    args_stat <- args_stat[1:(length(args_stat) - 1)]
    args_stat_i <- paste0("i_", args_stat)
    args_stat_rs <- paste0(args_stat, "_ls")

    if (base_fixed) par_base_ls <- h_ls <- list(NULL)

    for (i in 1:length(args_stat_rs)) {
        lenth_args_i <- length(eval(as.name(args_stat_rs[i])))

        if (lenth_args_i == 1) {
            assign(args_stat_i[i], rep(1L, n_regime))

        } else if (lenth_args_i == n_regime) {
            assign(args_stat_i[i], 1:n_regime)
        } else {
            stop("length of '", args_stat_rs[i], "' must be 1 or 'n_regime'")
        }
    }

    corr_ls <- vector("list", n_regime)

    for (i in 1:n_regime) {
        ind_i <- lapply(mget(args_stat_i), function(x)
            x[[i]])
        args_rs_i <- mget(args_stat_rs)
        names(args_rs_i) <- args_stat
        args_i <- Map(function(x, ind)
            x[[ind]], args_rs_i, ind_i)
        corr_ls[[i]] <-
            do.call(cor_stat, c(args_i, base_fixed = base_fixed))
    }
    return(corr_ls)
}
