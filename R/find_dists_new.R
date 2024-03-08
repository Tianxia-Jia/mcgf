#' Calculate (signed) distances between coordinates
#'
#' @param grid A matrix of 2D points, first column x/longitude, second column
#' y/latitude.
#' @param grid_new A matrix of 2D points, first column x/longitude, second column
#' y/latitude.
#' @param names_new Names of new locations.
#' @param longlat Logical, if TURE Great Circle (WGS84 ellipsoid) distance;
#' if FALSE, Euclidean distance.
#'
#' @keywords internal
#' @return List of signed distances between the new locations and the old grid.
.find_dists_new <- function(grid, grid_new, names = NULL, names_new = NULL,
                            longlat = TRUE) {
    grid_all <- rbind(grid, grid_new)
    n_var <- nrow(grid_all)
    names_all <- c(names, names_new)

    lat <- cbind(mean(grid[, 1]), grid_all[, 2])
    lat_dists <- sp::spDists(lat, longlat = longlat)
    rownames(lat_dists) <- colnames(lat_dists) <- names_all

    lon <- cbind(grid_all[, 1], mean(grid[, 2]))
    lon_dists <- sp::spDists(lon, longlat = longlat)
    rownames(lon_dists) <- colnames(lon_dists) <- names_all

    h <- sqrt(lon_dists^2 + lat_dists^2)
    rownames(h) <- colnames(h) <- names_all

    h1 <- matrix(0, ncol = n_var, nrow = n_var)
    for (i in 1:n_var) {
        for (j in 1:n_var) {
            h1[i, j] <- sign(grid_all[, 1][i] - grid_all[, 1][j]) *
                lon_dists[i, j]
        }
    }
    rownames(h1) <- colnames(h1) <- names_all

    h2 <- matrix(0, ncol = n_var, nrow = n_var)
    for (i in 1:n_var) {
        for (j in 1:n_var) {
            h2[i, j] <- sign(grid_all[, 2][i] - grid_all[, 2][j]) *
                lat_dists[i, j]
        }
    }
    rownames(h2) <- colnames(h2) <- names_all

    return(list(h = h, h1 = h1, h2 = h2))
}

#' Calculate (signed) distances between coordinates
#'
#' @inherit .find_dists_new return
#'
#' @param locations A matrix or data.frame of 2D points, the first column is
#' x/longitude, and the second column is y/latitude.
#' @param locations_new A matrix or data.frame of 2D points, the first column is
#' x/longitude, and the second column is y/latitude.
#' @param longlat Logical, if TURE Great Circle (WGS84 ellipsoid) distance;
#' if FALSE, Euclidean distance.
#'
#' @export
#'
#' @details
#' `locations` must be a matrix or data.frame containing 2 columns,
#' first column x/longitude, and second column y/latitude.The row names of
#' `locations` are used as the names of the locations.
#'
#' @examples
#' lon <- c(110, 120, 130)
#' lat <- c(50, 55, 60)
#' locations <- cbind(lon, lat)
#' rownames(locations) <- paste("Site", 1:3)
#' find_dists(locations)
#'
#' locations_new <- c(115, 55)
#' find_dists_new(locations, locations_new)
find_dists_new <- function(locations, locations_new, longlat = TRUE) {
    if (NCOL(locations) != 2) {
        stop("`locations` must contain 2 columns", call. = FALSE)
    }

    if (is.vector(locations_new)) {
        if (length(locations_new) != 2) {
            stop("incorrect dimension for `locations_new`", call. = FALSE)
        }
        locations_new <- matrix(locations_new, nrow = 1)
    }

    names <- rownames(locations)
    if (is.null(names)) {
        names <- seq_len(nrow(locations))
    }

    names_new <- rownames(locations_new)
    if (is.null(names_new)) {
        names_new <- paste0("New_", seq_len(nrow(locations_new)))
    }

    if (any(table(names_new) > 1)) {
        stop("duplicate row names found in `locations_new`", call. = FALSE)
    }

    if (any(names_new %in% names)) {
        stop("duplicate row names found in `locations_new` and `locations`",
            call. = FALSE
        )
    }

    dists_ls <- .find_dists_new(
        grid = locations, grid_new = locations_new,
        names = names, names_new = names_new,
        longlat = longlat
    )
    return(dists_ls)
}
