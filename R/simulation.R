#' Simulated Markov chain Gaussian field
#'
#' Simulated MCGF for 10 locations.
#' @name sim1
#' @format `sim1`: a list containing a data.frame with 1000 rows and 10 columns
#' and a list of distances
#'
#' @details
#' `sim1` contains a simulated MCGF for 10 locations. It is simulated with a
#' separable base model and a triangular Lagrangian model. The true parameters
#' for the base model are: \eqn{\text{nugget} = 0, c = 0.001, \gamma = 0.5,
#' a = 0.5, \alpha = 0.8}, and those for the Lagrangian model are:
#' \eqn{v1 = 200, v2 = 200, k = 2, \lambda = 0.2}
#'
#' @examples
#' # Code used to generate and fit `sim1` (not run)
#' # library(mcgf)
#' # set.seed(123)
#' # h <- rdists(10)
#' #
#' # N <- 1000
#' # lag <- 5
#' #
#' # par_base <- list(par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
#' #                  par_t = list(a = 0.5, alpha = 0.8))
#' # par_lagr <- list(v1 = 200, v2 = 200, k = 2)
#' #
#' # sim1 <- mcgf_sim(
#' #     N = N,
#' #     base = "sep",
#' #     lagrangian = "lagr_tri",
#' #     par_base = par_base,
#' #     par_lagr = par_lagr,
#' #     lambda = 0.2,
#' #     dists = h,
#' #     lag = lag
#' # )
#' # sim1 <- sim1[-c(1:(lag + 1)), ]
#' # rownames(sim1) <- 1:nrow(sim1)
#' #
#' # sim1 <- list(data = sim1, dists = h)
#' @family (simulated) datasets
"sim1"

#' Simulated regime-switching Markov chain Gaussian field
#'
#' Simulated RS-MCGF for 10 locations.
#' @name sim2
#' @format `sim2`: a list containing a data.frame with 1000 rows and 10 columns,
#' a list of distances, and a vector of regime labels.
#' @details
#' `sim2` contains a simulated RS-MCGF for 10 locations. It is simulated with
#' a regime-switching separable base model. The true parameters for the base
#' model are: \deqn{\text{Regime 1}: \text{nugget} = 0, c = 0.01, \gamma = 0.5,
#' a = 0.5, \alpha = 0.2,}
#' \deqn{\text{Regime 2}: \text{nugget} = 0, c = 0.04, \gamma = 0.5, a = 0.3,
#' \alpha = 0.9.}
#'
#' @examples
#' # Code used to generate and fit `sim2` (not run)
#' # library(mcgf)
#' # set.seed(123)
#' # h <- rdists(10)
#' #
#' # # simulate regimes
#' # K <- 2
#' # N <- 1000
#' # lag <- 5
#' #
#' # tran_mat <-
#' #     matrix(rnorm(K ^ 2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
#' # diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
#' # tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
#' # tran_mat
#' #
#' # regime <- rep(NA, N)
#' # regime[1] <- 1
#' #
#' # for (n in 2:(N)) {
#' #     regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
#' # }
#' # table(regime)
#' #
#' # # simulate RS MCGF
#' # par_base1 <- list(par_s = list(nugget = 0, c = 0.01, gamma = 0.5),
#' #                   par_t = list(a = 0.5, alpha = 0.2))
#' #
#' # par_base2 <- list(par_s = list(nugget = 0, c = 0.04, gamma = 0.5),
#' #                   par_t = list(a = 0.3, alpha = 0.9))
#' #
#' # sim2 <- mcgf_rs_sim(
#' #     N = N,
#' #     label = regime,
#' #     base_ls = list("sep"),
#' #     lagrangian_ls = list("none"),
#' #     par_base_ls = list(par_base1, par_base2),
#' #     lambda_ls = list(0.1, 0.3),
#' #     lag_ls = list(lag, lag),
#' #     dists_ls = list(h, h)
#' # )
#' # sim2 <- sim2[-c(1:(lag + 1)),]
#' # rownames(sim2) <- 1:nrow(sim2)
#' #
#' # sim2 <- list(data = sim2[,-1],
#' #              dists = h,
#' #              label = sim2[, 1])
#' @family (simulated) datasets
"sim2"

#' Simulated regime-switching Markov chain Gaussian field
#'
#' Simulated RS-MCGF for 20 locations.
#' @name sim3
#' @format `sim3`: a list containing a data.frame with 5000 rows and 20 columns
#' and a list of locations.
#'
#' @details
#' `sim3` contains a simulated RS-MCGF for 20 locations. It is simulated with
#' the same base model and a regime-switching Lagrangian model. The true
#' parameters for the base model are: \eqn{\text{nugget} = 0, c = 0.05,
#' \gamma = 0.5, a = 0.5, \alpha = 0.2}, and the true parameters for the
#' Lagrangian model are
#' \deqn{\text{Regime 1}: \lambda = 0.2, v_1 = -100, v_2 = 100, k = 2,}
#' \deqn{\text{Regime 1}: \lambda = 0.2, v_1 = 200, v_2 = 200, k = 2.}
#' For parameter estimation, the base model is assumed known and is used to
#' estimate the regime-switching prevailing winds.
#'
#' @examples
#' # Code used to generate and fit `sim3` (not run)
#' # library(mcgf)
#' # set.seed(123)
#' # h <- rdists(10)
#' #
#' # # simulate regimes
#' # K <- 2
#' # N <- 1000
#' # lag <- 5
#' #
#' # tran_mat <-
#' #     matrix(rnorm(K ^ 2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
#' # diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
#' # tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
#' # tran_mat
#' #
#' # regime <- rep(NA, N)
#' # regime[1] <- 1
#' #
#' # for (n in 2:(N)) {
#' #     regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1],])
#' # }
#' # table(regime)
#' #
#' # # simulate RS MCGF
#' # par_base <- list(par_s = list(nugget = 0, c = 0.05, gamma = 0.5),
#' #                  par_t = list(a = 0.5, alpha = 0.2))
#' #
#' # par_lagr1 <- list(v1 = -100, v2 = 100, k = 2)
#' # par_lagr2 <- list(v1 = 200, v2 = 200, k = 2)
#' #
#' # sim3 <- mcgf_rs_sim(
#' #     N = N,
#' #     label = regime,
#' #     base_ls = list("sep"),
#' #     lagrangian_ls = list("lagr_tri"),
#' #     par_base_ls = list(par_base),
#' #     par_lagr_ls = list(par_lagr1, par_lagr2),
#' #     lambda_ls = list(0.2, 0.2),
#' #     lag_ls = list(lag, lag),
#' #     dists_ls = list(h, h)
#' # )
#' # sim3 <- sim3[-c(1:(lag + 1)), ]
#' # rownames(sim3) <- 1:nrow(sim3)
#' #
#' # sim3 <- list(data = sim3[, -1],
#' #              dists = h,
#' #              label = sim3[, 1])
#' @family (simulated) datasets
"sim3"
