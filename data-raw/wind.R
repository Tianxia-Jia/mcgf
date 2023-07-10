## code to prepare `wind` dataset goes here
library(gstat)
data(wind)
time <- ISOdate(wind$year + 1900, wind$month, wind$day)
wind <- wind[, -c(1:3, 6)]
wind <- wind * 1852 / 3600
wind <- data.frame(time = time, wind)
wind <- wind[, c("time", "RPT", "VAL", "KIL", "SHA", "BIR", "DUB", "MUL", "CLA",
                 "CLO", "BEL", "MAL")]
usethis::use_data(wind, overwrite = TRUE)

library(sp)
lat <- as.numeric(char2dms(as.character(wind.loc[["Latitude"]])))
lon <- as.numeric(char2dms(as.character(wind.loc[["Longitude"]])))
wind_loc <- data.frame(lon = lon, lat = lat)
rownames(wind_loc) <- wind.loc$Code
wind_loc <- wind_loc[-12, ]
stat <- colnames(wind)[-1]
wind_loc <- wind_loc[match(stat, rownames(wind_loc)), ]
usethis::use_data(wind_loc, overwrite = TRUE)
