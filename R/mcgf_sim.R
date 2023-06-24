#' Simulate Markov chain Gaussian field
#'
#' @param N Sample size.
#' @param base Base model, `sep` or `fs` for now.
#' @param lagrangian Lagrangian model, "none" or `lagr_tri` for now.
#' @param par_base Parameters for the base model (symmetric).
#' @param par_lagr Parameters for the Lagrangian model.
#' @param Weight of the Lagrangian term, \eqn{\lambda\in[0, 1]}.
#' @param dists Distance matrices or arrays.
#' @param sd Standard Deviation for each location.
#' @param lag Time lag.
#' @param horizon Forecast horizon, default is 1.
#' @param init Initial samples, default is 0.
#' @param mu_c,mu_p Mean of current and past.
#' @param return_all Logical; if TRUE the joint covariance matrix, arrays of
#' distances and time lag are returned.
#'
#' @keywords internal
#'
#' @return Simulated Markov chain Gaussian field with user-specified covariance
#' structure. The simulation is done by kriging. The output data is in
#' space-wide format. `dists` must contain `h` for symmetric models, and `h1`
#' and `h2` for general stationary models. `horizon` controls forecasting
#' horizon. `sd`, `mu_c`, `init`, and `mu_p` must be vectors of appropriate
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
                      horizon,
                      init,
                      mu_c,
                      mu_p,
                      return_all = FALSE) {
    lag_max <- lag + horizon - 1

    n_var <- nrow(dists$h)
    n_block_row <- lag_max + 1

    stopifnot(N >= horizon + n_block_row)
    n_rounds <- ceiling((N - n_block_row) / horizon)

    u <- 0:lag_max
    dim_ar <- c(n_var, n_var, length(u))
    h_ar <- array(dists$h, dim = dim_ar)
    u_ar <- array(rep(u, each = n_var * n_var), dim = dim_ar)

    if (lagrangian == "lagr_tri") {
        h1_ar <- array(dists$h1, dim = dim_ar)
        h2_ar <- array(dists$h2, dim = dim_ar)

        cov_ar <-
            cor_stat(
                base = base,
                lagrangian = "lagr_tri",
                par_base = par_base,
                par_lagr = par_lagr,
                lambda = lambda,
                h = h_ar,
                h1 = h1_ar,
                h2 = h2_ar,
                u = u_ar,
                base_fixed = FALSE
            )
    } else {
        par_base <- append(par_base, list(h = h_ar, u = u_ar))
        cov_ar <- do.call(paste0("cor_", base), par_base)
    }

    for (i in 1:dim(cov_ar)[3]) {
        cov_ar[, , i] <- cor2cov(cov_ar[, , i], sd)
    }

    cov_mat_joint <- cov_joint(cov = cov_ar)

    ind_curr <- 1:(n_var * horizon)
    cov_mat_curr <- cov_mat_joint[ind_curr, ind_curr]
    cov_mat_curr_past <- cov_mat_joint[ind_curr, -ind_curr]
    cov_mat_past <- cov_mat_joint[-ind_curr, -ind_curr]
    lse <- cov_mat_curr_past %*% solve(cov_mat_past)

    X <- init

    for (n in 1:n_rounds) {
        X_past <- embed(tail(X, lag), lag)

        X_new_mean <- mu_c + lse %*% t(X_past - mu_p)
        X_new_cov <- cov_mat_curr - lse %*% t(cov_mat_curr_past)

        X_new <- mvnfast::rmvn(1, X_new_mean, X_new_cov)
        X_new <- matrix(X_new, ncol = n_var, byrow = T)
        X_new <- X_new[horizon:1, ]

        X <- rbind(X, X_new)
    }

    if (return_cov) {
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
#' @param N Sample size.
#' @param base Base model, `sep` or `fs` for now.
#' @param lagrangian Lagrangian model, "none" or `lagr_tri` for now.
#' @param par_base Parameters for the base model (symmetric).
#' @param par_lagr Parameters for the Lagrangian model.
#' @param Weight of the Lagrangian term, \eqn{\lambda\in[0, 1]}.
#' @param dists Distance matrices or arrays.
#' @param sd Standard Deviation for each location.
#' @param lag Time lag.
#' @param horizon Forecast horizon, default is 1.
#' @param init Initial samples, default is 0.
#' @param mu_c,mu_p Mean of current and past.
#' @param return_all Logical; if TRUE the joint covariance matrix, arrays of
#' distances and time lag are returned.
#'
#' @return Simulated Markov chain Gaussian field with user-specified covariance
#' structure. The simulation is done by kriging. The output data is in
#' space-wide format. `dists` must contain `h` for symmetric models and `h1`
#' and `h2` for general stationary models. `horizon` controls forecasting
#' horizon. `sd`, `mu_c`, `init`, and `mu_p` can be scalars.
#' @export
#'
#' @examples
#' 1
mcgf_sim <- function(N,
                     base = c("sep", "fs"),
                     lagrangian = c("none", "lagr_tri"),
                     par_base,
                     par_lagr,
                     lambda = 0,
                     dists,
                     sd = 1,
                     lag = 1,
                     horizon = 1,
                     init = 0,
                     mu_c = 0,
                     mu_p = 0,
                     return_all = FALSE) {

    if (N < horizon)
        stop("'N' must be no less than 'horizon'")

    if (is.null(dists$h))
        stop("missing 'h' is dists.")

    if (lagrangian == "lagr_tri") {

        if (is.null(dists$h1))
            stop("missing 'h1' is dists.")
        if (is.null(dists$h2))
            stop("missing 'h2' is dists.")
    }

    base <- match.arg(base)
    lagrangian <- match.arg(lagrangian)

    lag_max <- lag + horizon - 1
    n_var <- nrow(dist$h)
    n_block_row <- lag_max + 1

    if (length(sd) == 1) {
        sd <- rep(sd, n_var)
    } else {
        if (length(sd) == n_var)
            stop("length of 'sd' must be 1 or ", n_var)
    }

    if (length(init) == 1) {
        init <- matrix(init, nrow = n_block_row, ncol = n_var)
    } else {
        if (NROW(init) != n_block_row || NCOL(init) !=  n_var)
            stop("dim of 'n_var' must be 1 or ", n_block_row, " x ", n_var)
    }

    if (length(mu_c) == 1) {
        mu_c <- rep(mu_c, n_var * horizon)
    } else {
        if (length(mu_c) != n_var * horizon)
            stop("length of 'mu_c' must be 1 or ", n_var * horizon)
    }

    if (length(mu_p) == 1) {
        mu_p <- rep(mu_p, n_var * lag)
    } else {
        if (length(mu_p) != n_var * lag)
            stop("length of 'mu_p' must be 1 or ", n_var * lag)
    }

    n_rounds <- ceiling((N - n_block_row) / horizon)

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
