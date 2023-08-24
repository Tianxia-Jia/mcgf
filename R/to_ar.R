#' Convert to array
#'
#' @param h Distance matrix.
#' @param lag_max Maximum time lag.
#' @param u Logical; TRUE if u_ar needs to be outputted.
#'
#' @keywords internal
#' @return A list of arrays of h and u.
to_ar <- function(h, lag_max, u = TRUE) {
    dim_ar <- c(dim(h)[1:2], lag_max + 1)

    if (length(dim(h)) == 3) {
        if (any(dim_ar != dim(h)))
            stop('unmatching dimension of h.')
        h_ar <- h
    } else {
        h_ar <- array(h, dim = dim_ar)
    }

    if (u) {
        u_ar <- array(rep(0:lag_max, each = prod(dim(h)[1:2])), dim = dim_ar)
        return(list(h_ar = h_ar, u_ar = u_ar))
    } else {
        return(h_ar)
    }
}
