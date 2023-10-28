#' Print an `mcgf` object.
#'
#' @param x An `mcgf` object.
#' @param attr Attribute to be printed.
#' @param ... Optional arguments to print methods.
#'
#' @export
print.mcgf <- function(x, attr = ".Data", ...) {
    if (attr == ".Data") {
        print.data.frame(x, ...)
        cat(
            "\nOther attributes:",
            paste(names(attributes(x))[-c(1:3)], collapse = ", "),
            "\n"
        )
        return(invisible(NULL))
    } else {
        print(attr(x, attr, exact = TRUE))
        return(invisible(NULL))
    }
}
