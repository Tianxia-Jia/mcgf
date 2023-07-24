#' Ireland wind data, 1961-1978
#'
#' Daily average wind speeds for 1961-1978 at 11 synoptic meteorological
#' stations in the Republic of Ireland (Haslett and raftery 1989). Wind speeds
#' are in m/s. De-trended data sets are also provided.
#'
#' @name wind
#' @format `wind`: an object of class data.frame with 6574 rows and 12 columns.
#'
#' @details
#' The data were obtained from the **gstat** package, and were modified so that
#' the first column is the time stamps. Locations of the 11 stations are given
#' in `wind_loc`. `wind_train` and `wind_test` contain de-trended and
#' square-root transformed train (1961-1970) and test (1971-1978) data sets.
#' See Gneiting et al. (2006) for de-trending details. `wind_trend` contains
#' the estimated annual trend and station-wise mean from the training dataset.
#'
#' @references
#' Haslett, J. and Raftery, A. E. (1989). Space-time Modelling with Long-memory
#' Dependence: Assessing Ireland's Wind Power Resource (with Discussion).
#' Applied Statistics 38, 1-50.
#'
#' Gneiting, T., Genton, M., & Guttorp, P. (2006). Geostatistical Space-Time
#' Models, Stationarity, Separability, and Full Symmetry. In C&amp;H/CRC
#' Monographs on Statistics &amp; Applied Probability (pp. 151–175).
#' Chapman and Hall/CRC.
"wind"

#' Ireland wind data, 1961-1978
#' @rdname wind
#' @format `wind_loc`: an object of class data.frame with 12 rows and 2 columns.
"wind_loc"

#' #' Ireland wind data, 1961-1978
#' #' @rdname wind
#' #' @format `wind_train`: an object of class data.frame with 3650 rows and
#' #' 12 columns.
#' "wind_train"
#'
#' #' Ireland wind data, 1961-1978
#' #' @rdname wind
#' #' @format `wind_test`: an object of class data.frame with 2920 rows and
#' #' 12 columns.
#' "wind_test"
#'
#' #' Ireland wind data, 1961-1978
#' #' @rdname wind
#' #' @format `wind_trend`: an object of class list with first element of a
#' #' data.frame with 365 rows and 3 columns, and second element of a vector of
#' #' lenth 11.
#' "wind_trend"
