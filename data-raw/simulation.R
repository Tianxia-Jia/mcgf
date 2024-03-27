#------------------------------------------------------------------------------#
# Code for simulating a MCGF
#------------------------------------------------------------------------------#
library(mcgf)

set.seed(123)
x <- stats::rnorm(10, -110)
y <- stats::rnorm(10, 50)
locations <- cbind(x, y)
h <- find_dists(locations, longlat = TRUE)

N <- 1000
lag <- 5

par_base <- list(
    par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.8)
)
par_lagr <- list(v1 = 200, v2 = 200, k = 2)

sim1 <- mcgf_sim(
    N = N, base = "sep", lagrangian = "lagr_tri",
    par_base = par_base, par_lagr = par_lagr, lambda = 0.2,
    dists = h, lag = lag
)
sim1 <- sim1[-c(1:(lag + 1)), ]
rownames(sim1) <- 1:nrow(sim1)

sim1 <- list(data = sim1, locations = locations, dists = h)
usethis::use_data(sim1, overwrite = TRUE)

sim1_mcgf <- mcgf(sim1$data, locations = locations)
sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = lag)
sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = lag)

fit_spatial <- fit_base(sim1_mcgf,
    model = "spatial", lag = lag,
    par_init = c(c = 0.001, gamma = 0.5),
    par_fixed = c(nugget = 0)
)
fit_spatial$fit
# $par
#           c       gamma
# 0.001145201 0.500000000
#
# $objective
# [1] 1.385845
#
# $convergence
# [1] 0
#
# $iterations
# [1] 7
#
# $evaluations
# function gradient
#       18       17
#
# $message
# [1] "both X-convergence and relative convergence (5)"

fit_temporal <- fit_base(sim1_mcgf,
    model = "temporal", lag = lag,
    par_init = c(a = 0.3, alpha = 0.5)
)
fit_temporal$fit
# $par
#         a     alpha
# 0.6534102 0.7572530
#
# $objective
# [1] 0.004445817
#
# $convergence
# [1] 0
#
# $iterations
# [1] 16
#
# $evaluations
# function gradient
#       18       38
#
# $message
# [1] "relative convergence (4)"

sim1_mcgf <- add_base(sim1_mcgf,
    fit_s = fit_spatial, fit_t = fit_temporal,
    sep = TRUE
)

fit_sep <- fit_base(sim1_mcgf,
    model = "sep", lag = lag,
    par_init = c(c = 0.001, gamma = 0.5, a = 0.3, alpha = 0.5),
    par_fixed = c(nugget = 0)
)
fit_sep$fit
# $par
#           c       gamma           a       alpha
# 0.001139264 0.500000000 0.627518829 0.734384140
#
# $objective
# [1] 3.196846
#
# $convergence
# [1] 0
#
# $iterations
# [1] 25
#
# $evaluations
# function gradient
#       46      113
#
# $message
# [1] "relative convergence (4)"

sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep, old = TRUE)

fit_lagr <- fit_lagr(sim1_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 300, v2 = 300, lambda = 0.15),
    par_fixed = c(k = 2)
)
fit_lagr$fit
# $par
#      lambda          v1          v2
#   0.1625903 224.5776223 188.5158090
#
# $objective
# [1] 1.533448
#
# $convergence
# [1] 0
#
# $iterations
# [1] 34
#
# $evaluations
# function gradient
#       37      118
#
# $message
# [1] "relative convergence (4)"

sim1_mcgf_fitted <- add_lagr(sim1_mcgf, fit_lagr = fit_lagr)
model(sim1_mcgf_fitted, old = TRUE)
# ----------------------------------------
#                  Model
# ----------------------------------------
# - Time lag: 5
# - Scale of time lag: 1
# - Forecast horizon: 1
# ----------------------------------------
#             Old - not in use
# ----------------------------------------
# - Base-old model: sep
# - Parameters:
#           c       gamma      nugget           a       alpha
# 0.001145201 0.500000000 0.000000000 0.653410220 0.757253028
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls wls
#
# - Optimization function: nlminb nlminb
# ----------------------------------------
#                  Base
# ----------------------------------------
# - Base model: sep
# - Parameters:
#           c       gamma           a       alpha      nugget
# 0.001139264 0.500000000 0.627518829 0.734384140 0.000000000
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb
# ----------------------------------------
#               Lagrangian
# ----------------------------------------
# - Lagrangian model: lagr_tri
# - Parameters:
#      lambda          v1          v2           k
#   0.1625903 224.5776223 188.5158090   2.0000000
#
# - Fixed parameters:
# k
# 2
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb

#------------------------------------------------------------------------------#
# Code for simulating a RS MCGF in base
#------------------------------------------------------------------------------#
set.seed(123)
x <- stats::rnorm(10, -110)
y <- stats::rnorm(10, 50)
locations <- cbind(x, y)
h <- find_dists(locations, longlat = TRUE)

# simulate regimes
K <- 2
N <- 1000
lag <- 5

tran_mat <- matrix(rnorm(K^2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
tran_mat
#            [,1]       [,2]
# [1,] 0.94635675 0.05364325
# [2,] 0.06973429 0.93026571

regime <- rep(NA, N)
regime[1] <- 1

for (n in 2:(N)) {
    regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
}
table(regime)
# regime
#   1   2
# 567 433

# simulate RS MCGF
par_base1 <- list(
    par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)

par_base2 <- list(
    par_s = list(nugget = 0, c = 0.004, gamma = 0.5),
    par_t = list(a = 0.3, alpha = 0.9)
)

sim2 <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"), lagrangian_ls = list("none"),
    par_base_ls = list(par_base1, par_base2),
    lambda_ls = list(0.1, 0.3),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
sim2 <- sim2[-c(1:(lag + 1)), ]
rownames(sim2) <- 1:nrow(sim2)

sim2 <- list(
    data = sim2[, -1], locations = locations, dists = h,
    label = sim2[, 1]
)
usethis::use_data(sim2, overwrite = TRUE)

# Parameter estimation
sim2_mcgf <- mcgf_rs(sim2$data, locations = locations, label = sim2$label)
sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = lag)
sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = lag)

fit_spatial <- fit_base(sim2_mcgf,
    lag_ls = lag, model_ls = "spatial",
    par_init_ls = list(c(c = 0.00005, gamma = 0.5)),
    par_fixed_ls = list(c(nugget = 0))
)
lapply(fit_spatial[1:2], function(x) x$fit)
# $`Regime 1`
# $`Regime 1`$par
#           c       gamma
# 0.001496696 0.500000000
#
# $`Regime 1`$objective
# [1] 3.325991
#
# $`Regime 1`$convergence
# [1] 0
#
# $`Regime 1`$iterations
# [1] 60
#
# $`Regime 1`$evaluations
# function gradient
#       81      144
#
# $`Regime 1`$message
# [1] "both X-convergence and relative convergence (5)"
#
#
# $`Regime 2`
# $`Regime 2`$par
#           c       gamma
# 0.004074408 0.500000000
#
# $`Regime 2`$objective
# [1] 1.106692
#
# $`Regime 2`$convergence
# [1] 0
#
# $`Regime 2`$iterations
# [1] 48
#
# $`Regime 2`$evaluations
# function gradient
#       61      118
#
# $`Regime 2`$message
# [1] "both X-convergence and relative convergence (5)"

fit_temporal <- fit_base(sim2_mcgf,
    lag_ls = lag, model_ls = "temporal",
    par_init_ls = list(
        list(a = 0.8, alpha = 0.8),
        list(a = 0.1, alpha = 0.1)
    )
)
lapply(fit_temporal[1:2], function(x) x$fit)
# $`Regime 1`
# $`Regime 1`$par
#         a     alpha
# 0.4375028 0.3020300
#
# $`Regime 1`$objective
# [1] 0.007620233
#
# $`Regime 1`$convergence
# [1] 0
#
# $`Regime 1`$iterations
# [1] 14
#
# $`Regime 1`$evaluations
# function gradient
#       20       35
#
# $`Regime 1`$message
# [1] "relative convergence (4)"
#
#
# $`Regime 2`
# $`Regime 2`$par
#         a     alpha
# 0.3309062 0.8745089
#
# $`Regime 2`$objective
# [1] 0.005785562
#
# $`Regime 2`$convergence
# [1] 0
#
# $`Regime 2`$iterations
# [1] 22
#
# $`Regime 2`$evaluations
# function gradient
#       27       52
#
# $`Regime 2`$message
# [1] "both X-convergence and relative convergence (5)"

sim2_mcgf <- add_base(sim2_mcgf,
    fit_s_ls = fit_spatial,
    fit_t_ls = fit_temporal, sep = TRUE
)
model(sim2_mcgf, model = "base")
# ----------------------------------------
#                  Model
# ----------------------------------------
# - Time lag: 5, 5
# - Scale of time lag: 1
# - Forecast horizon: 1
# ----------------------------------------
#                  Base
# ----------------------------------------
# - Regime switching: spatial: TRUE, temporal: TRUE
# --------------------
#       Regime 1
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma      nugget           a       alpha
# 0.001496696 0.500000000 0.000000000 0.437502762 0.302029968
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls wls
#
# - Optimization function: nlminb nlminb
# --------------------
#       Regime 2
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma      nugget           a       alpha
# 0.004074408 0.500000000 0.000000000 0.330906238 0.874508851
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls wls
#
# - Optimization function: nlminb nlminb

fit_sep <- fit_base(sim2_mcgf,
    lag_ls = lag, model_ls = "sep",
    par_init_ls = list(list(
        c = 0.00005, gamma = 0.5,
        a = 0.5, alpha = 0.5
    )),
    par_fixed_ls = list(c(nugget = 0))
)
lapply(fit_sep[1:2], function(x) x$fit)
# $`Regime 1`
# $`Regime 1`$par
#           c       gamma           a       alpha
# 0.001539941 0.500000000 0.452841529 0.345950120
#
# $`Regime 1`$objective
# [1] 7.62732
#
# $`Regime 1`$convergence
# [1] 0
#
# $`Regime 1`$iterations
# [1] 62
#
# $`Regime 1`$evaluations
# function gradient
#       90      290
#
# $`Regime 1`$message
# [1] "relative convergence (4)"
#
#
# $`Regime 2`
# $`Regime 2`$par
#           c       gamma           a       alpha
# 0.004107168 0.500000000 0.328388504 0.853896071
#
# $`Regime 2`$objective
# [1] 3.302677
#
# $`Regime 2`$convergence
# [1] 0
#
# $`Regime 2`$iterations
# [1] 60
#
# $`Regime 2`$evaluations
# function gradient
#       80      279
#
# $`Regime 2`$message
# [1] "relative convergence (4)"

sim2_mcgf_fitted <- add_base(sim2_mcgf, fit_base_ls = fit_sep, old = TRUE)
model(sim2_mcgf_fitted, model = "base", old = TRUE)
# ----------------------------------------
#                  Model
# ----------------------------------------
# - Time lag: 5, 5
# - Scale of time lag: 1
# - Forecast horizon: 1
# ----------------------------------------
#             Old - not in use
# ----------------------------------------
# - Regime switching: spatial: TRUE, temporal: TRUE
# --------------------
#       Regime 1
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma      nugget           a       alpha
# 0.001496696 0.500000000 0.000000000 0.437502762 0.302029968
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls wls
#
# - Optimization function: nlminb nlminb
# --------------------
#       Regime 2
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma      nugget           a       alpha
# 0.004074408 0.500000000 0.000000000 0.330906238 0.874508851
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls wls
#
# - Optimization function: nlminb nlminb
# ----------------------------------------
#                  Base
# ----------------------------------------
# - Regime switching: TRUE
# --------------------
#       Regime 1
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma           a       alpha      nugget
# 0.001539941 0.500000000 0.452841529 0.345950120 0.000000000
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb
# --------------------
#       Regime 2
# --------------------
# - Base model: sep
# - Parameters:
#           c       gamma           a       alpha      nugget
# 0.004107168 0.500000000 0.328388504 0.853896071 0.000000000
#
# - Fixed parameters:
# nugget
#      0
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb
#------------------------------------------------------------------------------#
# Code for simulating a RS MCGF in lagr
#------------------------------------------------------------------------------#
set.seed(123)
x <- stats::rnorm(10, -110)
y <- stats::rnorm(10, 50)
locations <- cbind(x, y)
h <- find_dists(locations, longlat = TRUE)

# simulate regimes
K <- 2
N <- 1000
lag <- 5

tran_mat <- matrix(rnorm(K^2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
tran_mat
# [,1]       [,2]
# [1,] 0.94635675 0.05364325
# [2,] 0.06973429 0.93026571

regime <- rep(NA, N)
regime[1] <- 1

for (n in 2:(N)) {
    regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
}
table(regime)
# regime
#   1   2
# 567 433

# simulate RS MCGF
par_base <- list(
    par_s = list(nugget = 0, c = 0.05, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)

par_lagr1 <- list(v1 = -100, v2 = 100, k = 2)
par_lagr2 <- list(v1 = 200, v2 = 200, k = 2)

sim3 <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"), lagrangian_ls = list("lagr_tri"),
    par_base_ls = list(par_base),
    par_lagr_ls = list(par_lagr1, par_lagr2),
    lambda_ls = list(0.2, 0.2),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
sim3 <- sim3[-c(1:(lag + 1)), ]
rownames(sim3) <- 1:nrow(sim3)

sim3 <- list(
    data = sim3[, -1], locations = locations, dists = h,
    label = sim3[, 1]
)
usethis::use_data(sim3, overwrite = TRUE)

# Parameter estimation

sim3_mcgf <- mcgf_rs(sim3$data, locations = locations, label = sim3$label)
sim3_mcgf <- add_acfs(sim3_mcgf, lag_max = lag)
sim3_mcgf <- add_ccfs(sim3_mcgf, lag_max = lag)

fit_fs <- fit_base(sim3_mcgf,
    lag_ls = lag, model_ls = "fs", rs = FALSE,
    par_init_ls = list(list(beta = 0)),
    par_fixed_ls = list(list(
        nugget = 0, c = 0.05, gamma = 0.5,
        a = 0.5, alpha = 0.2
    ))
)
fit_fs[[1]]$fit$par <- 0

sim3_mcgf <- add_base(sim3_mcgf, fit_base_ls = fit_fs)
model(sim3_mcgf, model = "base")
# ----------------------------------------
#                  Model
# ----------------------------------------
# - Time lag: 5
# - Scale of time lag: 1
# - Forecast horizon: 1
# ----------------------------------------
#                  Base
# ----------------------------------------
# - Regime switching: FALSE
# - Base model: fs
# - Parameters:
#   beta nugget      c  gamma      a  alpha
#   0.00   0.00   0.05   0.50   0.50   0.20
#
# - Fixed parameters:
# nugget      c  gamma      a  alpha
#   0.00   0.05   0.50   0.50   0.20
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb

fit_lagr_rs <- fit_lagr(sim3_mcgf,
    model_ls = list("lagr_tri"),
    par_init_ls = list(
        list(v1 = -50, v2 = 50),
        list(v1 = 100, v2 = 100)
    ),
    par_fixed_ls = list(list(lambda = 0.2, k = 2))
)
lapply(fit_lagr_rs[1:2], function(x) x$fit)
# $`Regime 1`
# $`Regime 1`$par
#        v1        v2
# -106.5181  119.4434
#
# $`Regime 1`$objective
# [1] 5.137749
#
# $`Regime 1`$convergence
# [1] 0
#
# $`Regime 1`$iterations
# [1] 26
#
# $`Regime 1`$evaluations
# function gradient
#       29       61
#
# $`Regime 1`$message
# [1] "relative convergence (4)"
#
#
# $`Regime 2`
# $`Regime 2`$par
#       v1       v2
# 210.1156 242.3082
#
# $`Regime 2`$objective
# [1] 5.902847
#
# $`Regime 2`$convergence
# [1] 0
#
# $`Regime 2`$iterations
# [1] 29
#
# $`Regime 2`$evaluations
# function gradient
#       31       77
#
# $`Regime 2`$message
# [1] "relative convergence (4)"

sim3_mcgf_fitted <- add_lagr(sim3_mcgf, fit_lagr_ls = fit_lagr_rs)
model(sim3_mcgf_fitted)
# ----------------------------------------
#                  Model
# ----------------------------------------
# - Time lag: 5, 5
# - Scale of time lag: 1
# - Forecast horizon: 1
# ----------------------------------------
#                  Base
# ----------------------------------------
# - Regime switching: FALSE
# - Base model: fs
# - Parameters:
#   beta nugget      c  gamma      a  alpha
#   0.00   0.00   0.05   0.50   0.50   0.20
#
# - Fixed parameters:
# nugget      c  gamma      a  alpha
#   0.00   0.05   0.50   0.50   0.20
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb
# ----------------------------------------
#               Lagrangian
# ----------------------------------------
# - Regime switching: TRUE
# --------------------
#       Regime 1
# --------------------
# - Lagrangian model: lagr_tri
# - Parameters:
#        v1        v2    lambda         k
# -106.5181  119.4434    0.2000    2.0000
#
# - Fixed parameters:
# lambda      k
#    0.2    2.0
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb
# --------------------
#       Regime 2
# --------------------
# - Lagrangian model: lagr_tri
# - Parameters:
#       v1       v2   lambda        k
# 210.1156 242.3082   0.2000   2.0000
#
# - Fixed parameters:
# lambda      k
#    0.2    2.0
#
# - Parameter estimation method: wls
#
# - Optimization function: nlminb

fit_lagr_rs_mle <- fit_lagr(sim3_mcgf,
    model_ls = list("lagr_tri"),
    par_init_ls = list(
        list(v1 = -50, v2 = 50),
        list(v1 = 100, v2 = 100)
    ),
    method = "mle",
    par_fixed_ls = list(list(lambda = 0.2, k = 2))
)
lapply(fit_lagr_rs_mle[1:2], function(x) x$fit)
# $`Regime 1`
# $`Regime 1`$par
#        v1        v2
# -106.1742  101.4039
#
# $`Regime 1`$objective
# [1] 1144.1
#
# $`Regime 1`$convergence
# [1] 0
#
# $`Regime 1`$iterations
# [1] 12
#
# $`Regime 1`$evaluations
# function gradient
#       15       32
#
# $`Regime 1`$message
# [1] "relative convergence (4)"
#
#
# $`Regime 2`
# $`Regime 2`$par
#       v1       v2
# 215.6930 189.3384
#
# $`Regime 2`$objective
# [1] 973.9486
#
# $`Regime 2`$convergence
# [1] 0
#
# $`Regime 2`$iterations
# [1] 20
#
# $`Regime 2`$evaluations
# function gradient
#       25       48
#
# $`Regime 2`$message
# [1] "relative convergence (4)"
