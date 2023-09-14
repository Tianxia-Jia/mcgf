#' Generic functions for calculating joint covariance/correlation matrix for mcgf
#' objects
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Joint correlation/covariance matrix.
#' @export
#' @family {functions related to cor and cov}
ccov <- function(x, ...) {
    UseMethod("ccov")
}

#' Covariance/correlation for joint distribution of an `mcgf` object
#'
#' @param x An `mcgf` object.
#' @param model Which model to use. One of `all`, `base`, and `empirical`.
#' @param cor Logical; if TRUE, correlation is outputted.
#' @param ... Additional arguments. Not in use.
#'
#' @return Joint vovariance/correlation matrix.
#' @export
ccov.mcgf <- function(x, model = c("all", "base", "empirical"),
                     cor = FALSE, ...) {

    model <- match.arg(model)

    horizon <- attr(x, "horizon", exact = TRUE)
    lag <- attr(x, "lag", exact = TRUE)
    n_var <- ncol(dists(x)$h)
    empirical <- FALSE

    if (model == "empirical") {
        cor_ar <- ccfs(x)[, , 1:(horizon + lag)]
        empirical <- TRUE

    } else if (model == "base") {
        cor_ar <- attr(x, "base_res", exact = TRUE)$cor_base
    } else {
        cor_ar <- attr(x, "lagr_res", exact = TRUE)$cor_lagr
    }

    if (!cor) {
        cor_ar <- cor2cov_ar(cor_ar, sds(x), empirical = empirical)
    }

    cov_mat <- cov_joint(cor_ar)
    loc_name <- colnames(dists(x)$h)
    name_hori <- paste0(loc_name,
                        rep(paste0(":t+", (horizon - 1):0), each = n_var))
    name_lag <- paste0(loc_name,
                       rep(paste0(":lag", 1:(lag)), each = n_var))
    colnames(cov_mat) <- rownames(cov_mat) <- c(name_hori, name_lag)

    return(cov_mat)
}

#' Covariance/correlation for joint distribution of an `mcgf_rs`object
#'
#' @param x An `mcgf` object.
#' @param model Which model to use. One of `all`, `base`, and `empirical`.
#' @param cor Logical; if TRUE, correlation is returned
#' @param ... Additional arguments. Not in use.
#'
#' @return A list of joint vovariance/correlation matrix.
#' @export
ccov.mcgf_rs <- function(x, model = c("all", "base", "empirical"),
                         cor = FALSE, ...) {

    model <- match.arg(model)
    horizon <- attr(x, "horizon", exact = TRUE)
    lag_ls <- attr(x, "lag", exact = TRUE)
    n_var <- ncol(dists(x)$h)

    empirical <- FALSE
    lvs <- levels((attr(x, "label", exact = TRUE)))
    n_regime <- length(lvs)

    if (model == "empirical") {
        cor_ar_ls <- Map(function(x, y) x[, , 1:(horizon + y)],
                         ccfs(x)$ccfs_rs, lag_ls)
        empirical <- TRUE
    } else if (model == "base") {
        cor_ar_ls <- lapply(attr(x, "base_res", exact = TRUE),
                            function(x) x$cor_base)
    } else {
        cor_ar_ls <- lapply(attr(x, "lagr_res", exact = TRUE),
                            function(x) x$cor_lagr)
    }

    if (!cor) {
        cor_ar_ls <- Map(cor2cov_ar, cor_ar_ls, sds(x)$sds_rs,
                         empirical = empirical)
    }

    cov_mat_ls <- lapply(cor_ar_ls, cov_joint)
    loc_name <- colnames(dists(x)$h)
    name_hori <- paste0(loc_name, rep(paste0(":t+", (horizon -
                                                         1):0), each = n_var))
    for (i in 1:n_regime) {
        name_lag <- paste0(loc_name, rep(paste0(":lag", 1:lag_ls[[i]]),
                                         each = n_var))
        colnames(cov_mat_ls[[i]]) <- rownames(cov_mat_ls[[i]]) <-
            c(name_hori, name_lag)
    }
    return(cov_mat_ls)
}
