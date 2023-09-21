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
    "time", "VAL", "BEL", "CLA", "SHA", "RPT", "BIR", "MUL", "MAL",
    "KIL", "CLO", "DUB"
)]
colnames(wind)[colnames(wind) == "RPT"] <- "ROC"
usethis::use_data(wind, overwrite = TRUE)

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
usethis::use_data(wind_loc, overwrite = TRUE)

#------------------------------------------------------------------------------#
# Code for preparing wind_de datasets
#------------------------------------------------------------------------------#
# library(tidyr)
# library(dplyr)
# library(tibble)
#
# # Remove leap day and take square-root
# is_leap <- month(wind$time) == 2 & day(wind$time) == 29
# wind_de <- wind[!ind_leap, ]
# wind_de[, -1] <- sqrt(wind_de[, -1])
#
# # Split into training and test
# is_train <- year(wind_de$time) <= 1970
# wind_train <- wind_de[is_train, ]
# wind_test <- wind_de[!is_train, ]
#
# # Estimate annual trend
# avg_train <- tibble(month = month(wind_train$time),
#                     day = day(wind_train$time),
#                     ws = rowMeans(wind_train[,-1]))
# trend_train <- avg_train %>%
#     group_by(month, day) %>%
#     summarise(trend = mean(ws))
#
# # Subtract annual trend
# trend_train2 <- left_join(avg_train, trend_train, by = c("month", "day"))$trend
# wind_train[, -1] <- wind_train[, -1] - trend_train2
#
# # Subtract station-wise mean
# mean_train <- colMeans(wind_train[, -1])
# wind_train[, -1] <- t(apply(wind_train[, -1], 1,
#                                function(x) x - mean_train))
# wind_trend <- list(annual = as.data.frame(trend_train),
#                    mean = mean_train)
#
# # Do the same for the test data
# avg_test <- tibble(month = month(wind_test$time),
#                    day = day(wind_test$time))
# trend_train3 <- left_join(avg_test, trend_train, by = c("month", "day"))$trend
# wind_test[, -1] <- wind_test[, -1] - trend_train3
# wind_test[, -1] <- t(apply(wind_test[, -1], 1,
#                               function(x) x - mean_train))
# usethis::use_data(wind_train, overwrite = TRUE)
# usethis::use_data(wind_test, overwrite = TRUE)
# usethis::use_data(wind_trend, overwrite = TRUE)
