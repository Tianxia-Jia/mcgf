#' Generic function for displaying fitted models for `mcgf` objects
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Details of the fitted models.
#' @export
#'
#' @details
#' Refer to [`model.mcgf()`] and [`model.mcgf_rs()`] for more details.
model <- function(x, ...) {
    UseMethod("model")
}

#' Display fitted models for an `mcgf` or `mcgf_rs` object
#'
#' @name model.mcgf
#' @param x An mcgf object.
#' @param model Which model to display.
#' @param old Logical; TRUE if the old model needs to be printed.
#' @param print_model Logical; TRUE if time lag and forecast horizon need to be
#' printed.
#' @param ... Additional arguments. Not in use.
#'
#' @return None (invisible `NULL`).
#' @export
#'
#' @details
#' For `mcgf` and `mcgf_rs` objects, [`model()`] displays the fitted models and
#' their parameters. When `old = TRUE`, the old model is printed as well. Note
#' that the old model is not used for parameter estimation or for kriging.
#'
#' @family functions on fitting a mcgf
#' @family functions on fitting a mcgf_rs
model.mcgf <- function(x, model = c("all", "base", "lagrangian"), old = FALSE,
                       print_model = TRUE, ...) {
    model <- match.arg(model)

    if (print_model) {
        cat("----------------------------------------\n")
        cat("                 Model\n")
        cat("----------------------------------------\n")
        cat("- Time lag:", attr(x, "lag", exact = TRUE), "\n")
        cat("- Scale of time lag:", attr(x, "scale_time", exact = TRUE), "\n")
        cat(
            "- Forecast horizon:",
            attr(x, "horizon", exact = TRUE),
            "\n"
        )
    }

    if (old) {
        base_res <- attr(x, "base_res_old", exact = TRUE)

        cat("----------------------------------------\n")
        cat("            Old - not in use\n")
        cat("----------------------------------------\n")
        cat("- Base-old model:", attr(x, "base_old", exact = TRUE), "\n")
        cat("- Parameters:\n")
        print(unlist(base_res$par_base))
        cat("\n- Fixed parameters:\n")
        print(unlist(base_res$par_fixed))
        cat("\n- Parameter estimation method:", base_res$method_base, "\n")
        cat("\n- Optimization function:", base_res$optim_fn, "\n")
    }

    if (model == "base") {
        base_res <- attr(x, "base_res", exact = TRUE)

        cat("----------------------------------------\n")
        cat("                 Base\n")
        cat("----------------------------------------\n")
        cat("- Base model:", attr(x, "base", exact = TRUE), "\n")
        cat("- Parameters:\n")
        print(unlist(base_res$par_base))
        cat("\n- Fixed parameters:\n")
        print(unlist(base_res$par_fixed))
        cat("\n- Parameter estimation method:", base_res$method_base, "\n")
        cat("\n- Optimization function:", base_res$optim_fn, "\n")
    } else if (model == "lagrangian") {
        lagr_res <- attr(x, "lagr_res", exact = TRUE)

        cat("----------------------------------------\n")
        cat("              Lagrangian\n")
        cat("----------------------------------------\n")
        cat("- Lagrangian model:", attr(x, "lagr", exact = TRUE), "\n")
        cat("- Parameters:\n")
        print(unlist(lagr_res$par_lagr))
        cat("\n- Fixed parameters:\n")
        print(unlist(lagr_res$par_fixed))
        cat("\n- Parameter estimation method:", lagr_res$method_lagr, "\n")
        cat("\n- Optimization function:", lagr_res$optim_fn, "\n")
    } else {
        model.mcgf(x, "base", print_model = FALSE)
        model.mcgf(x, "lagrangian", print_model = FALSE)
    }
    invisible(NULL)
}

#' @rdname model.mcgf
#' @export
model.mcgf_rs <- function(x, model = c("all", "base", "lagrangian"),
                          old = FALSE, print_model = TRUE, ...) {
    model <- match.arg(model)

    if (print_model) {
        cat("----------------------------------------\n")
        cat("                 Model\n")
        cat("----------------------------------------\n")
        cat(
            "- Time lag:",
            paste(attr(x, "lag", exact = TRUE), collapse = ", "), "\n"
        )
        cat("- Scale of time lag:", attr(x, "scale_time", exact = TRUE), "\n")
        cat("- Forecast horizon:", attr(x, "horizon", exact = TRUE), "\n")
    }

    if (old) {
        base_rs <- attr(x, "base_rs_old", exact = TRUE)
        base_res <- attr(x, "base_res_old", exact = TRUE)

        cat("----------------------------------------\n")
        cat("            Old - not in use\n")
        cat("----------------------------------------\n")

        if (length(base_rs) == 2 && sum(base_rs) == 0) {
            cat("- Regime switching: FALSE\n")
            cat("- Base model:", attr(x, "base_old", exact = TRUE), "\n")

            cat("- Parameters:\n")
            print(unlist(base_res$par_base))
            cat("\n- Fixed parameters:\n")
            print(unlist(base_res$par_fixed))
            cat(
                "\n- Parameter estimation method:",
                base_res$method_base, "\n"
            )
            cat("\n- Optimization function:", base_res$optim_fn, "\n")
        } else if (length(base_rs) == 1 && !base_rs) {
            cat("- Regime switching: FALSE\n")
            cat("- Base model:", attr(x, "base_old", exact = TRUE), "\n")

            cat("- Parameters:\n")
            print(unlist(base_res$par_base))
            cat("\n- Fixed parameters:\n")
            print(unlist(base_res$par_fixed))
            cat(
                "\n- Parameter estimation method:",
                base_res$method_base, "\n"
            )
            cat("\n- Optimization function:", base_res$optim_fn, "\n")
        } else {
            if (length(base_rs) == 1) {
                cat("- Regime switching:", base_rs, "\n")
            } else {
                cat(
                    "- Regime switching:",
                    paste0(names(base_rs), ": ", base_rs, collapse = ", "),
                    "\n"
                )
            }

            lvs <- levels(attr(x, "label", exact = TRUE))
            n_regime <- length(lvs)

            for (i in 1:n_regime) {
                r_name <- paste0("Regime ", lvs[i])
                cat("--------------------\n")
                cat("     ", r_name, "\n")
                cat("--------------------\n")
                cat(
                    "- Base model:",
                    attr(x, "base_old", exact = TRUE)[[i]], "\n"
                )
                cat("- Parameters:\n")
                print(unlist(base_res[[i]]$par_base))
                cat("\n- Fixed parameters:\n")
                print(unlist(base_res[[i]]$par_fixed))
                cat(
                    "\n- Parameter estimation method:",
                    base_res[[i]]$method_base, "\n"
                )
                cat(
                    "\n- Optimization function:",
                    base_res[[i]]$optim_fn, "\n"
                )
            }
        }
    }

    if (model == "base") {
        base_rs <- attr(x, "base_rs", exact = TRUE)
        base_res <- attr(x, "base_res", exact = TRUE)

        cat("----------------------------------------\n")
        cat("                 Base\n")
        cat("----------------------------------------\n")

        if (length(base_rs) == 2 && sum(base_rs) == 0) {
            cat("- Regime switching: FALSE\n")
            cat("- Base model:", attr(x, "base", exact = TRUE), "\n")

            cat("- Parameters:\n")
            print(unlist(base_res$par_base))
            cat("\n- Fixed parameters:\n")
            print(unlist(base_res$par_fixed))
            cat("\n- Parameter estimation method:", base_res$method_base, "\n")
            cat("\n- Optimization function:", base_res$optim_fn, "\n")
        } else if (length(base_rs) == 1 && !base_rs) {
            cat("- Regime switching: FALSE\n")
            cat("- Base model:", attr(x, "base", exact = TRUE), "\n")

            cat("- Parameters:\n")
            print(unlist(base_res$par_base))
            cat("\n- Fixed parameters:\n")
            print(unlist(base_res$par_fixed))
            cat("\n- Parameter estimation method:", base_res$method_base, "\n")
            cat("\n- Optimization function:", base_res$optim_fn, "\n")
        } else {
            if (length(base_rs) == 1) {
                cat("- Regime switching:", base_rs, "\n")
            } else {
                cat(
                    "- Regime switching:",
                    paste0(names(base_rs), ": ", base_rs, collapse = ", "),
                    "\n"
                )
            }

            lvs <- levels(attr(x, "label", exact = TRUE))
            n_regime <- length(lvs)

            for (i in 1:n_regime) {
                r_name <- paste0("Regime ", lvs[i])
                cat("--------------------\n")
                cat("     ", r_name, "\n")
                cat("--------------------\n")
                cat("- Base model:", attr(x, "base", exact = TRUE)[[i]], "\n")
                cat("- Parameters:\n")
                print(unlist(base_res[[i]]$par_base))
                cat("\n- Fixed parameters:\n")
                print(unlist(base_res[[i]]$par_fixed))
                cat(
                    "\n- Parameter estimation method:",
                    base_res[[i]]$method_base, "\n"
                )
                cat("\n- Optimization function:", base_res[[i]]$optim_fn, "\n")
            }
        }
    } else if (model == "lagrangian") {
        lagr_rs <- attr(x, "lagr_rs", exact = TRUE)
        lagr_res <- attr(x, "lagr_res", exact = TRUE)

        cat("----------------------------------------\n")
        cat("              Lagrangian\n")
        cat("----------------------------------------\n")

        if (is.null(lagr_rs)) {
            cat("- Lagrangian model:", attr(x, "lagr", exact = TRUE), "\n")
            cat("- Parameters:\n")
            print(unlist(lagr_res$par_lagr))
            cat("\n- Fixed parameters:\n")
            print(unlist(lagr_res$par_fixed))
            cat("\n- Parameter estimation method:", lagr_res$method_lagr, "\n")
            cat("\n- Optimization function:", lagr_res$optim_fn, "\n")
        } else if (!lagr_rs) {
            cat("- Regime switching: FALSE\n")
            cat("- Lagrangian model:", attr(x, "lagr", exact = TRUE), "\n")
            cat("- Parameters:\n")
            print(unlist(lagr_res$par_lagr))
            cat("\n- Fixed parameters:\n")
            print(unlist(lagr_res$par_fixed))
            cat("\n- Parameter estimation method:", lagr_res$method_lagr, "\n")
            cat("\n- Optimization function:", lagr_res$optim_fn, "\n")
        } else {
            cat("- Regime switching:", lagr_rs, "\n")

            lvs <- levels(attr(x, "label", exact = TRUE))
            n_regime <- length(lvs)

            for (i in 1:n_regime) {
                r_name <- paste0("Regime ", lvs[i])
                cat("--------------------\n")
                cat("     ", r_name, "\n")
                cat("--------------------\n")
                cat(
                    "- Lagrangian model:",
                    attr(x, "lagr", exact = TRUE)[[i]], "\n"
                )
                cat("- Parameters:\n")
                print(unlist(lagr_res[[i]]$par_lagr))
                cat("\n- Fixed parameters:\n")
                print(unlist(lagr_res[[i]]$par_fixed))
                cat(
                    "\n- Parameter estimation method:",
                    lagr_res[[i]]$method_lagr, "\n"
                )
                cat("\n- Optimization function:", lagr_res[[i]]$optim_fn, "\n")
            }
        }
    } else {
        model.mcgf_rs(x, "base", print_model = FALSE)
        model.mcgf_rs(x, "lagrangian", print_model = FALSE)
    }
    return(invisible(NULL))
}
