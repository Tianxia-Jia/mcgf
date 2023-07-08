#' Print an `mcgf` object.
#'
#' @param x An `mcgf` object.
#' @param ... Optional arguments to print methods.
#' @param attr Attribute to be printed.
#'
#' @export
#' @family {functions related to the class}
print.mcgf <- function(x, attr = ".Data", ...) {

    if (attr == ".Data") {
        print.data.frame(x, ...)
        cat("\nOther attributes:",
            paste(names(attributes(x))[-c(1:3)], collapse = ", "), "\n")
        return(invisible(NULL))
    }

    if (attr == "acfs") {
        print(attr(x, "acfs"))
        return(invisible(NULL))
    }
}
