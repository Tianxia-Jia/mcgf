#' Generic function for standard deviations for each column
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Standard deviations.
#' @export
#' @family {functions related to the class}
sds <- function(x, ...) {
    UseMethod("sds")
}

#' @rdname ccfs.mcgf
#' @param value Cross-correlations.
sds.mcgf <- function(x) {
    return(attr(x, "sds"))
}

#' @rdname add_ccfs.mcgf
#' @param value A vector of standard deviations for all stations.
`sds<-` <- function(x, value) {

    if (!is.vector(x, mode = "numeric"))
        stop("`value` must be a numeric vector.")

    if (length(value) != ncol(x))
        stop("length of `value` must be the same as number of columns.")

    attr(x, "sds") <- value
}
