#' Generic function for displaying fitted models for mcgf objects
#'
#' @param x An **R** object.
#' @param ... Additional parameters or attributes.
#'
#' @return Fitted models for an mcgf object.
#' @export
#' @family {functions related to model fitting}
model <- function(x, ...) {
    UseMethod("model")
}

#' Display fitted models for an mcgf object
#'
#' @param x An mcgf object.
#' @param model Which model to display.
#' @param old Logical; TRUE if the old model needs to be printed.
#' @param print_model Logical; TRUE if time lag and forecast horizon need to be
#' printed.
#' @param ... Additional arguments. Not in use.
#'
#' @return None (invisible `NULL`).
#' @export
#' @family {functions related to model fitting}
model.mcgf <- function(x, model = c("all", "base", "lagrangian"), old = FALSE,
                       print_model = TRUE, ...) {

        model <- match.arg(model)

        if (print_model) {
            cat('----------------------------------------\n')
            cat('                 Model\n')
            cat('----------------------------------------\n')
            cat("- Time lag:", attr(x, "lag", exact = TRUE), "\n")
            cat("- Forecast horizon:",
                attr(x, "horizon", exact = TRUE),
                "\n")
        }

        if (model == "base") {
            cat('----------------------------------------\n')
            cat('                 Base\n')
            cat('----------------------------------------\n')
            cat("- Base model:", attr(x, "base", exact = TRUE), "\n")

            cat("- Parameters:\n")
            print(unlist(attr(x, "base_res", exact = TRUE)$par_base))
            cat("\n- Fixed parameters:\n")
            print(unlist(attr(x, "base_res", exact = TRUE)$par_fixed))
            cat(
                "\n- Parameter estimation method:",
                attr(x, "base_res", exact = TRUE)$method_base,
                "\n"
            )
            cat("\n- Optimization function:",
                attr(x, "base_res", exact = TRUE)$optim_fn,
                "\n")
        } else if (model == "lagrangian") {
            cat('----------------------------------------\n')
            cat('              Lagrangian\n')
            cat('----------------------------------------\n')
            cat("- Lagrangian model:", attr(x, "lagr", exact = TRUE), "\n")
            cat("- Parameters:\n")
            print(unlist(attr(x, "lagr_res", exact = TRUE)$par_lagr))
            cat("\n- Fixed parameters:\n")
            print(unlist(attr(x, "lagr_res", exact = TRUE)$par_fixed))
            cat(
                "\n- Parameter estimation method:",
                attr(x, "lagr_res", exact = TRUE)$method_lagr,
                "\n"
            )
            cat("\n- Optimization function:",
                attr(x, "lagr_res", exact = TRUE)$optim_fn,
                "\n")
        } else {
            if (old) {
                cat('----------------------------------------\n')
                cat('            Old - not in use\n')
                cat('----------------------------------------\n')
                cat("- Base-old model:",
                    attr(x, "base_old", exact = TRUE),
                    "\n")
                cat("- Parameters:\n")
                print(unlist(attr(x, "base_res_old", exact = TRUE)$par_base))
                cat("\n- Fixed parameters:\n")
                print(unlist(attr(x, "base_res_old", exact = TRUE)$par_fixed))
                cat(
                    "\n- Parameter estimation method:",
                    attr(x, "base_res_old", exact = TRUE)$method_base,
                    "\n"
                )
                cat(
                    "\n- Optimization function:",
                    attr(x, "base_res_old", exact = TRUE)$optim_fn,
                    "\n"
                )
            }

            model.mcgf(x, "base", print_model = FALSE)
            model.mcgf(x, "lagrangian", print_model = FALSE)
        }
        invisible(NULL)
}
