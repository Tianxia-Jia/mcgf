#------------------------------------------------------------------------------#
# Code for preparing wind dataset
#------------------------------------------------------------------------------#
library(gstat)
library(lubridate)

data(wind)
time <- as.Date(paste(wind$year + 1900, wind$month, wind$day, sep = "-"))
wind <- wind[, -c(1:3, 6)]
wind <- wind * 1852 / 3600
wind <- data.frame(time = time, wind)
wind <- wind[, c(
    "time",
    "VAL",
    "BEL",
    "CLA",
    "SHA",
    "RPT",
    "BIR",
    "MUL",
    "MAL",
    "KIL",
    "CLO",
    "DUB"
)]
colnames(wind)[colnames(wind) == "RPT"] <- "ROC"

#------------------------------------------------------------------------------#
# Code for preparing wind_loc dataset
#------------------------------------------------------------------------------#
library(sp)

wind.loc$Code <- as.character(wind.loc$Code)
wind.loc$Code[wind.loc$Code == "RPT"] <- "ROC"
lat <- as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
lon <- as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
wind_loc <- data.frame(lon = lon, lat = lat)
rownames(wind_loc) <- wind.loc$Code
wind_loc <- wind_loc[-12, ]
stat <- colnames(wind)[-1]
wind_loc <- wind_loc[match(stat, rownames(wind_loc)), ]

wind <- list(data = wind, locations = wind_loc)
usethis::use_data(wind, overwrite = TRUE)

#------------------------------------------------------------------------------#
# Code for preparing wind_de datasets
#------------------------------------------------------------------------------#
library(tidyr)
library(dplyr)
library(tibble)

# Remove leap day and take square-root
ind_leap <- month(wind$data$time) == 2 & day(wind$data$time) == 29
wind_de <- wind$data[!ind_leap, ]
wind_de[, -1] <- sqrt(wind_de[, -1])

# Split into training and test
is_train <- year(wind_de$time) <= 1970
wind_train <- wind_de[is_train, ]
wind_test <- wind_de[!is_train, ]

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
trend_train2 <-
    left_join(avg_train, trend_train, by = c("month", "day"))$trend
wind_train[, -1] <- wind_train[, -1] - trend_train2

# Subtract station-wise mean
mean_train <- colMeans(wind_train[, -1])
wind_train[, -1] <- t(apply(
    wind_train[, -1], 1,
    function(x) {
        x - mean_train
    }
))
wind_trend <- list(
    annual = as.data.frame(trend_train),
    mean = mean_train
)

# Do the same for the test data
avg_test <- tibble(
    month = month(wind_test$time),
    day = day(wind_test$time)
)
trend_train3 <-
    left_join(avg_test, trend_train, by = c("month", "day"))$trend
wind_test[, -1] <- wind_test[, -1] - trend_train3
wind_test[, -1] <- t(apply(
    wind_test[, -1], 1,
    function(x) x - mean_train
))

#------------------------------------------------------------------------------#
# Code for fitting covariance model on the de-trended wind dataset
#------------------------------------------------------------------------------#
wind_mcgf <- mcgf(wind_train[, -1], locations = wind$locations)
wind_mcgf <- add_acfs(wind_mcgf, lag_max = 3)
wind_mcgf <- add_ccfs(wind_mcgf, lag_max = 3)

fit_spatial <- fit_base(
    wind_mcgf,
    model = "spatial",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001),
    par_fixed = c(gamma = 0.5)
)
fit_spatial$fit
# $par
#           c      nugget
# 0.001326192 0.049136042
#
# $objective
# [1] 1.742162
#
# $convergence
# [1] 0
#
# $iterations
# [1] 9
#
# $evaluations
# function gradient
#       23       24
#
# $message
# [1] "relative convergence (4)"

fit_temporal <- fit_base(
    wind_mcgf,
    model = "temporal",
    lag = 3,
    par_init = c(a = 0.5, alpha = 0.5)
)
fit_temporal$fit
# $par
#         a     alpha
# 0.9774157 0.8053363
#
# $objective
# [1] 0.0006466064
#
# $convergence
# [1] 0
#
# $iterations
# [1] 12
#
# $evaluations
# function gradient
#       15       31
#
# $message
# [1] "both X-convergence and relative convergence (5)"

wind_mcgf <-
    add_base(wind_mcgf,
        fit_s = fit_spatial,
        fit_t = fit_temporal,
        sep = T
    )

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
# $par
#      beta
# 0.6232458
#
# $objective
# [1] 2.858384
#
# $convergence
# [1] 0
#
# $iterations
# [1] 5
#
# $evaluations
# function gradient
#        6        7
#
# $message
# [1] "relative convergence (4)"

wind_mcgf <- add_base(wind_mcgf, fit_base = fit_fs, old = T)
model(wind_mcgf, old = T)

fit_lagr <- fit_lagr(wind_mcgf,
    model = "lagr_tri",
    par_init = c(v2 = 200, lambda = 0.1),
    par_fixed = c(v1 = 0, k = 2)
)
fit_lagr$fit
# $par
#       lambda           v2
#   0.03535575 329.03455269
#
# $objective
# [1] 2.798634
#
# $convergence
# [1] 0
#
# $iterations
# [1] 17
#
# $evaluations
# function gradient
#       20       48
#
# $message
# [1] "relative convergence (4)"

fit_lagr2 <- fit_lagr(wind_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 100, v2 = 100, lambda = 0.1),
    par_fixed = c(k = 2)
)
fit_lagr2$fit
# $par
#       lambda           v1           v2
#   0.03919917 202.74574080 236.28860761
#
# $objective
# [1] 2.779051
#
# $convergence
# [1] 0
#
# $iterations
# [1] 28
#
# $evaluations
# function gradient
#       32      145
#
# $message
# [1] "relative convergence (4)"

fit_lagr3 <- fit_lagr(wind_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 100, v2 = 100, lambda = 0.1, k = 2)
)
fit_lagr3$fit
# $par
#       lambda           v1           v2            k
#   0.04000252 194.56772556 255.26984386   1.94074140
#
# $objective
# [1] 2.778418
#
# $convergence
# [1] 0
#
# $iterations
# [1] 96
#
# $evaluations
# function gradient
#      141      648
#
# $message
# [1] "relative convergence (4)"
