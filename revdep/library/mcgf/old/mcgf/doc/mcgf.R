## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>"
)

## ----load data----------------------------------------------------------------
library(mcgf)
data(wind)
head(wind$data)
head(wind$locations)

## ----leap, message = F--------------------------------------------------------
# install.packages("lubridate")
library(lubridate)
ind_leap <- month(wind$data$time) == 2 & day(wind$data$time) == 29
wind_de <- wind$data[!ind_leap, ]
wind_de[, -1] <- sqrt(wind_de[, -1])

## ----split--------------------------------------------------------------------
is_train <- year(wind_de$time) <= 1970
wind_train <- wind_de[is_train, ]
wind_test <- wind_de[!is_train, ]

## ----annual, message = F------------------------------------------------------
# install.packages("dplyr")
library(dplyr)
# Estimate annual trend
avg_train <- tibble(
    month = month(wind_train$time),
    day = day(wind_train$time),
    ws = rowMeans(wind_train[, -1])
)
trend_train <- avg_train %>%
    group_by(month, day) %>%
    summarise(trend = mean(ws))

# Subtract annual trend
trend_train2 <- left_join(avg_train, trend_train, by = c("month", "day"))$trend
wind_train[, -1] <- wind_train[, -1] - trend_train2

## ----station------------------------------------------------------------------
# Subtract station-wise mean
mean_train <- colMeans(wind_train[, -1])
wind_train[, -1] <- sweep(wind_train[, -1], 2, mean_train)
wind_trend <- list(
    annual = as.data.frame(trend_train),
    mean = mean_train
)

## ----test---------------------------------------------------------------------
avg_test <- tibble(
    month = month(wind_test$time),
    day = day(wind_test$time)
)
trend_train3 <-
    left_join(avg_test, trend_train, by = c("month", "day"))$trend
wind_test[, -1] <- wind_test[, -1] - trend_train3
wind_test[, -1] <- sweep(wind_test[, -1], 2, mean_train)

## ----mcgf---------------------------------------------------------------------
wind_mcgf <- mcgf(wind_train[, -1], locations = wind$locations)
wind_mcgf <- add_acfs(wind_mcgf, lag_max = 3)
wind_mcgf <- add_ccfs(wind_mcgf, lag_max = 3)

## ----acfs---------------------------------------------------------------------
acfs(wind_mcgf)

## ----ccfs---------------------------------------------------------------------
ccfs(wind_mcgf)

## ----spatial, message=F-------------------------------------------------------
fit_spatial <- fit_base(
    wind_mcgf,
    model = "spatial",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001),
    par_fixed = c(gamma = 0.5)
)
fit_spatial$fit

## -----------------------------------------------------------------------------
library(Rsolnp)
solnp2 <- function(pars, fun, lower, upper, ...) {
    solnp(pars, fun, LB = lower, UB = upper, ...)
}
fit_spatial2 <- fit_base(
    wind_mcgf,
    model = "spatial",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001),
    par_fixed = c(gamma = 0.5),
    optim_fn = "other",
    other_optim_fn = "solnp2"
)
fit_spatial2$fit

## ----temporal-----------------------------------------------------------------
fit_temporal <- fit_base(
    wind_mcgf,
    model = "temporal",
    lag = 3,
    par_init = c(a = 0.5, alpha = 0.5)
)
fit_temporal$fit

## ----sep----------------------------------------------------------------------
wind_mcgf <- add_base(wind_mcgf,
    fit_s = fit_spatial,
    fit_t = fit_temporal,
    sep = T
)

## -----------------------------------------------------------------------------
fit_sep <- fit_base(
    wind_mcgf,
    model = "sep",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001, a = 0.5, alpha = 0.5),
    par_fixed = c(gamma = 0.5)
)
fit_sep$fit

## ----fs-----------------------------------------------------------------------
par_sep <- c(fit_spatial$fit$par, fit_temporal$fit$par, gamma = 0.5)
fit_fs <-
    fit_base(
        wind_mcgf,
        model = "fs",
        lag = 3,
        par_init = c(beta = 0),
        par_fixed = par_sep
    )
fit_fs$fit

## ----base---------------------------------------------------------------------
wind_mcgf <- add_base(wind_mcgf, fit_base = fit_fs, old = TRUE)
model(wind_mcgf, model = "base", old = TRUE)

## ----westerly-----------------------------------------------------------------
fit_lagr <- fit_lagr(wind_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 200, lambda = 0.1),
    par_fixed = c(v2 = 0, k = 2)
)
fit_lagr$fit

## ----stat---------------------------------------------------------------------
wind_mcgf <- add_lagr(wind_mcgf, fit_lagr = fit_lagr)
model(wind_mcgf, old = TRUE)

## ----krige--------------------------------------------------------------------
krige_emp <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "empirical",
    interval = TRUE
)
krige_base <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "base",
    interval = TRUE
)
krige_stat <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "all",
    interval = TRUE
)

## ----rmse---------------------------------------------------------------------
# RMSE
rmse_emp <- sqrt(colMeans((wind_test[, -1] - krige_emp$fit)^2, na.rm = T))
rmse_base <- sqrt(colMeans((wind_test[, -1] - krige_base$fit)^2, na.rm = T))
rmse_stat <- sqrt(colMeans((wind_test[, -1] - krige_stat$fit)^2, na.rm = T))
rmse <- c(
    "Empirical" = mean(rmse_emp),
    "Base" = mean(rmse_base),
    "STAT" = mean(rmse_emp)
)
rmse

## ----MAE----------------------------------------------------------------------
mae_emp <- colMeans(abs(wind_test[, -1] - krige_emp$fit), na.rm = T)
mae_base <- colMeans(abs(wind_test[, -1] - krige_base$fit), na.rm = T)
mae_stat <- colMeans(abs(wind_test[, -1] - krige_stat$fit), na.rm = T)
mae <- c(
    "Empirical" = mean(mae_emp),
    "Base" = mean(mae_base),
    "STAT" = mean(mae_stat)
)
mae

## ----popi---------------------------------------------------------------------
# POPI
popi_emp <- colMeans(
    wind_test[, -1] < krige_emp$lower | wind_test[, -1] > krige_emp$upper,
    na.rm = T
)
popi_base <- colMeans(
    wind_test[, -1] < krige_base$lower | wind_test[, -1] > krige_base$upper,
    na.rm = T
)
popi_stat <- colMeans(
    wind_test[, -1] < krige_stat$lower | wind_test[, -1] > krige_stat$upper,
    na.rm = T
)
popi <- c(
    "Empirical" = mean(popi_emp),
    "Base" = mean(popi_base),
    "STAT" = mean(popi_stat)
)
popi

