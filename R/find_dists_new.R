#' Calculate (signed) distances between coordinates
#'
#' @param grid A matrix of 2D points, first column x/longitude, second column
#' y/latitude.
#' @param grid_new A matrix of 2D points, first column x/longitude, second column
#' y/latitude.
#' @param names Names of locations.
#' @param names_new Names of new locations.
#' @param longlat Logical, if TURE Great Circle (WGS84 ellipsoid) distance;
#' if FALSE, Euclidean distance.
#' @param origin Optional; used when `longlat` is TRUE. An integer index
#' indicating the reference location which will be used as the origin.
#' @param return_grid Logical; used when `longlat` is TRUE. If TRUE the mapped
#' coordinates on a 2D plane is returned.
#' @param lon_ref Reference longitude when computing the longitudinal distances.
#' Default is the mean of longitudes in `grid`.
#' @param lat_ref Reference latitude when computing the latitudinal distances.
#' Default is the mean of latitudes in `grid`.
#'
#' @keywords internal
#' @return List of signed distances between the new locations and the old grid.
.find_dists_new <- function(grid, grid_new, names = NULL, names_new = NULL,
                            longlat = TRUE, origin = 1L, return_grid = FALSE,
                            lon_ref = NULL, lat_ref = NULL) {
    grid_all <- rbind(grid, grid_new)
    n_var <- nrow(grid_all)
    names_all <- c(names, names_new)

    lon_ref <- ifelse(is.null(lon_ref), mean(grid[, 1]), lon_ref)
    lat_ref <- ifelse(is.null(lat_ref), mean(grid[, 2]), lat_ref)

    lon <- cbind(grid_all[, 1], lat_ref)
    lon_dists <- sp::spDists(lon, longlat = longlat)
    rownames(lon_dists) <- colnames(lon_dists) <- names_all

    lat <- cbind(lon_ref, grid_all[, 2])
    lat_dists <- sp::spDists(lat, longlat = longlat)
    rownames(lat_dists) <- colnames(lat_dists) <- names_all

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

    dists <- list(h = h, h1 = h1, h2 = h2)

    if (longlat) {
        grid_2d <- cbind(h1[, origin], h2[, origin])
        dists <- .find_dists(grid_2d, names = names_all, longlat = FALSE)

        if (return_grid) {
            return(list(dists = dists, grid = grid_2d, origin = origin))
        } else {
            return(dists)
        }
    } else {
        return(dists)
    }
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
#' @param origin Optional; used when `longlat` is TRUE. An integer index
#' indicating the reference location from `locations` which will be used as
#' the origin. Same `origin` from `find_dists` must be used to ensure
#' consistancy between outputs from `find_dists` and `find_dists_new`.
#' @param return_grid Logical; used when `longlat` is TRUE. If TRUE the mapped
#' coordinates on a 2D plane for all locations is returned.
#' @param ... Optional arguments passed to [`.find_dists_new()`].
#'
#' @return A list of distance matrices for all locations. If `return_grid` is
#' TRUE, a list consists of a list of distance matrices for all locations,
#' the mapped 2D grid for all locations, and the origin is returned.
#'
#' @export
#'
#' @details
#' `locations` and `locations_new` must be matrices or data.frames containing
#' 2 columns, first column x/longitude, and second column y/latitude. The row
#' names of `locations` and `locations_new` are used as the names of the
#' locations.
#'
#' If `longlat` is TRUE, the original coordinates are mapped to a 2D Euclidean
#' plane given the reference location from `locations`. First, the Great Circle
#' (WGS84 ellipsoid) signed distance matrices are calculated, where the original
#' latitudes are replaced by the the mean of latitudes in `locations` to find
#' the signed longitudinal distances and the original longitudes are replaced by
#' the the mean of longitudes in `locations` to find the signed latitudinal
#' distances. Then given the index of a reference location `origin`, a new set
#' of coordinates in a 2D plane is generated where the coordinates are
#' determined by the signed distances between the locations and the reference
#' location. Finally distance matrices of the new coordinates for all stations
#' are outputted.
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
find_dists_new <- function(locations, locations_new, longlat = TRUE,
                           origin = 1L, return_grid = FALSE, ...) {
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

    locations_new <- as.data.frame(locations_new)
    names(locations_new)[1:2] <- colnames(locations)[1:2]

    dists_ls <- .find_dists_new(
        grid = locations, grid_new = locations_new,
        names = names, names_new = names_new,
        longlat = longlat, origin = origin,
        return_grid = return_grid, ...
    )
    return(dists_ls)
}
