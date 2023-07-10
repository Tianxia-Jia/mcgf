#' Convert to array
#'
#' @param x Distance matrix.
#' @param lag_max Maximum time lag.
#'
#' @keywords internal
#' @return A list of arrays of h and u.
to_ar <- function(h, lag_max) {
    dim_ar <- c(dim(h), lag_max + 1)
    h_ar <- array(h, dim = dim_ar)
    u_ar <- array(rep(0:lag_max, each = prod(dim(h))), dim = dim_ar)
    return(list(h_ar = h_ar, u_ar = u_ar))
}
