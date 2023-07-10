#' Ireland wind data, 1961-1978
#'
#' Daily average wind speeds for 1961-1978 at 11 synoptic meteorological
#' stations in the Republic of Ireland (Haslett and raftery 1989). Wind speeds
#' are in m/s.
#'
#' @name wind
#' @aliases `wind_loc`
#' @format `wind`: an object of class data.frame with 6574 rows and 12 columns.
#'
#' @details
#' The data were obtained from the **gstat** package, and were modified so that
#' the first column is the time stamps. Locations of the 11 stations are given
#' in `wind_loc`.
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
