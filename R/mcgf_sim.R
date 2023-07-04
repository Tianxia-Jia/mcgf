#' Simulate Markov chain Gaussian field
#'
#' @param N Sample size.
#' @param base Base model, `sep` or `fs` for now.
#' @param lagrangian Lagrangian model, "none" or `lagr_tri` for now.
#' @param par_base Parameters for the base model (symmetric).
#' @param par_lagr Parameters for the Lagrangian model.
#' @param lambda Weight of the Lagrangian term, \eqn{\lambda\in[0, 1]}.
#' @param dists Distance matrices or arrays.
#' @param sd Standard deviation for each location.
#' @param lag Time lag.
#' @param horizon Forecast horizon, default is 1.
#' @param init Initial samples, default is 0.
#' @param mu_c,mu_p Means of current and past.
#' @param return_all Logical; if TRUE the joint covariance matrix, arrays of
#' distances and time lag are returned.
#'
#' @keywords internal
#'
#' @return Simulated Markov chain Gaussian field with user-specified covariance
#' structure. The simulation is done by kriging. The output data is in
#' space-wide format. `dists` must contain `h` for symmetric models, and `h1`
#' and `h2` for general stationary models. `horizon` controls forecasting
#' horizon. `sd`, `mu_c`, `mu_p`, and `init` must be vectors of appropriate
#' sizes.
.mcgf_sim <- function(N,
                      base,
                      lagrangian,
                      par_base,
                      par_lagr,
                      lambda,
                      dists,
                      sd,
                      lag,
                      horizon = 1,
                      init = 0,
                      mu_c,
                      mu_p,
                      return_all = FALSE) {

    lag_max <- lag + horizon - 1
    n_var <- nrow(dists$h)
    n_block_row <- lag_max + 1

    n_rounds <- ceiling((N - n_block_row) / horizon)

    u <- 0:lag_max
    dim_ar <- c(n_var, n_var, length(u))
    h_ar <- array(dists$h, dim = dim_ar)
    u_ar <- array(rep(u, each = n_var * n_var), dim = dim_ar)

    if (lagrangian == "none") {
        cov_ar <- cor_stat(
            base = base,
            lagrangian = lagrangian,
            par_base = par_base,
            lambda = lambda,
            h = h_ar,
            u = u_ar,
            base_fixed = FALSE
        )
    } else {
        h1_ar <- array(dists$h1, dim = dim_ar)
        h2_ar <- array(dists$h2, dim = dim_ar)

        cov_ar <- cor_stat(
            base = base,
            lagrangian = lagrangian,
            par_base = par_base,
            par_lagr = par_lagr,
            lambda = lambda,
            h = h_ar,
            h1 = h1_ar,
            h2 = h2_ar,
            u = u_ar,
            base_fixed = FALSE
        )
    }

    for (i in 1:dim(cov_ar)[3]) {
        cov_ar[, , i] <- cor2cov(cov_ar[, , i], sd)
    }

    cov_mat_joint <- cov_joint(cov = cov_ar)

    ind_curr <- 1:(n_var * horizon)
    cov_mat_curr <- cov_mat_joint[ind_curr, ind_curr]
    cov_mat_curr_past <- cov_mat_joint[ind_curr, -ind_curr]
    cov_mat_past <- cov_mat_joint[-ind_curr, -ind_curr]
    LSE <- cov_mat_curr_past %*% solve(cov_mat_past)
    X_new_cov <- cov_mat_curr - LSE %*% t(cov_mat_curr_past)

    X <- init

    for (n in 1:n_rounds) {
        X_past <- embed(tail(X, lag), lag)
        X_new_mean <- mu_c + LSE %*% t(X_past - mu_p)

        X_new <- mvnfast::rmvn(1, X_new_mean, X_new_cov)
        X_new <- matrix(X_new, ncol = n_var, byrow = T)
        X_new <- X_new[horizon:1, ]
        X <- rbind(X, X_new)
    }

    rownames(X) <- 1:nrow(X)
    colnames(X) <- colnames(dists$h)

    if (return_all) {
        if (lagrangian == "lagr_tri") {
            par <- list(
                cov_mat = cov_mat_joint,
                dists = list(
                    h = h_ar,
                    h1 = h1_ar,
                    h2 = h2_ar
                ),
                u = u_ar
            )
            return(list(X = X, par = par))
        }

    } else {
        return(X = X)
    }
}

#' Simulate Markov chain Gaussian field
#'
#' @inherit .mcgf_sim params details return
#'
#' @export
#' @examples
#' par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
#' par_t <- list(a = 1, alpha = 0.5)
#' par_base <- list(par_s = par_s, par_t = par_t)
#' par_lagr <- list(v1 = 5, v2 = 10)
#' h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
#' h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
#' h <- sqrt(h1 ^ 2 + h2 ^ 2)
#' dists <- list(h = h, h1 = h1, h2 = h2)
#'
#' set.seed(123)
#' X <- mcgf_sim(N = 1000, base = "sep", lagrangian = "lagr_tri", lambda = 0.5,
#'               par_base = par_base, par_lagr = par_lagr, dists = dists)
#' # plot.ts(X)
#'
#' @family {simulations of Markov chain Gaussian fields}
mcgf_sim <- function(N,
                     base = c("sep", "fs"),
                     lagrangian = c("none", "lagr_tri", "lagr_askey"),
                     par_base,
                     par_lagr,
                     lambda,
                     dists,
                     sd = 1,
                     lag = 1,
                     horizon = 1,
                     init = 0,
                     mu_c = 0,
                     mu_p = 0,
                     return_all = FALSE) {

    if (N < horizon)
        stop("'N' must be no less than 'horizon'", call. = FALSE)

    if (is.null(dists$h))
        stop("missing 'h' in dists.", call. = FALSE)

    if (lagrangian != "none") {

        if (is.null(dists$h1))
            stop("missing 'h1' in dists.", call. = FALSE)
        if (is.null(dists$h2))
            stop("missing 'h2' in dists.", call. = FALSE)
    }

    lag_max <- lag + horizon - 1
    n_var <- nrow(dists$h)
    n_block_row <- lag_max + 1

    if (N < horizon + n_block_row) {
        warning("'N' must be no less than ", horizon + n_block_row)
        N <- horizon + n_block_row
    }

    if (length(init) == 1) {
        init <- matrix(init, nrow = n_block_row, ncol = n_var)
    } else {
        if (NROW(init) != n_block_row || NCOL(init) !=  n_var)
            stop("dim of 'n_var' must be 1 or ", n_block_row, " x ", n_var, ".",
                 call. = FALSE)
    }

    sd <- check_length(x = sd, length = n_var, name = "sd")
    mu_c <- check_length(x = mu_c, length = n_var * horizon, name = "mu_c")
    mu_p <- check_length(x = mu_p, length = n_var * lag, name = "mu_p")

    res <- .mcgf_sim(
        N = N,
        base = base,
        lagrangian = lagrangian,
        par_base = par_base,
        par_lagr = par_lagr,
        lambda = lambda,
        dists = dists,
        sd = sd,
        lag = lag,
        horizon = horizon,
        init = init,
        mu_c = mu_c,
        mu_p = mu_p,
        return_all = return_all
    )
    return(res)
}
