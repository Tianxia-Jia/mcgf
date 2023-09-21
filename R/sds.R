#' Generic function for standard deviations for each column
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @details
#' Refer to [`sds.mcgf()`] for more details.
#'
#' @export
sds <- function(x, ...) {
    UseMethod("sds")
}

#' Extract, calculate, or assign standard deviations for an `mcgf` or
#' `mcgf_rs` object.
#'
#' @name sds.mcgf
#'
#' @param x An `mcgf` or `mcgf_rs` object.
#' @param ... Additional parameters or attributes. Not in use.
#'
#' @return [`sds()`] returns empirical (regime-switching) standard deviations.
#'
#' @details
#' For `mcgf` objects, [`sds()`] extracts or computes the empirical standard
#' deviations. The output is a vector of sds.
#'
#' For `mcgf_rs` objects, [`sds()`] extracts or computes the regime-switching
#' empirical standard deviations. The output is a list of vectors of sds. Each
#' element in the list corresponds to the sds for a regime.
#'
#' [`sds<-`] assigns `sds` to `x`. Use [`add_ccfs()`] to add both `ccfs` and
#' `sds` to `x`.
#'
#' @export
sds.mcgf <- function(x, ...) {
    sds <- attr(x, "sds", exact = TRUE)

    if (!is.null(sds)) {
        return(sds)
    } else {
        sds <- apply(x, 2, stats::sd)
        return(sds)
    }
}

#' @rdname sds.mcgf
#'
#' @param replace Logical; if TRUE, `sds` are recalculated.
#' @export
sds.mcgf_rs <- function(x, replace = FALSE, ...) {
    if (replace) {
        label <- attr(x, "label", exact = TRUE)
        sds_rs <- sd_rs(x = x, label = label)
        sds <- sds.mcgf(x)
        return(list(sds = sds, sds_rs = sds_rs))
    } else {
        sds <- attr(x, "sds", exact = TRUE)
        if (!is.null(sds)) {
            return(sds)
        } else {
            return(sds.mcgf_rs(x, replace = TRUE))
        }
    }
}

#' Calculate standard deviation for each location under each regime.
#'
#' @param x A `data.frame` or `matrix`.
#' @param label A vector of regime labels. Its length must be the same as
#' the number rows in `x`.
#'
#' @return A list of standard deviations for each regime.
#' @export
sd_rs <- function(x, label) {
    if (!is.factor(label)) {
        label <- as.factor(label)
    }

    lvs <- levels(label)

    n_reg <- length(lvs)
    sd_ls <- lapply(1:n_reg, function(i) {
        apply(x[label == lvs[i], ], 2, stats::sd)
    })
    names(sd_ls) <- paste0("Regime ", levels(label))

    return(sd_ls)
}

#' @rdname sds.mcgf
#' @param value A vector (or list of vectors) of standard deviations for all
#' stations (under each regime and combined).
#' @export
`sds<-` <- function(x, value) {
    attr(x, "sds") <- value
    return(x)
}
