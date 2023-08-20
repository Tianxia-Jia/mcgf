#' Convert to array
#'
#' @param h Distance matrix.
#' @param lag_max Maximum time lag.
#' @param u Logical; TRUE if u_ar needs to be outputted.
#'
#' @keywords internal
#' @return A list of arrays of h and u.
to_ar <- function(h, lag_max, u = TRUE) {
    dim_ar <- c(dim(h), lag_max + 1)
    h_ar <- array(h, dim = dim_ar)

    if (u) {
        u_ar <- array(rep(0:lag_max, each = prod(dim(h))), dim = dim_ar)
        return(list(h_ar = h_ar, u_ar = u_ar))
    } else {
        return(h_ar)
    }
}
