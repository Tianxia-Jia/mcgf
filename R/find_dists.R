#' Calculate (signed) distances between coordinates
#'
#' @param grid A matrix of 2D points, first column x/longitude, second column
#' y/latitude.
#' @param names Names of locations.
#' @param longlat Logical, if TURE Great Circle (WGS84 ellipsoid) distance;
#' if FALSE, Euclidean distance.
#'
#' @keywords internal
#' @return List of signed distances.
.find_dists <- function(grid, names = NULL, longlat = TRUE) {
    n_var <- nrow(grid)

    lat <- cbind(mean(grid[, 1]), grid[, 2])
    lat_dists <- sp::spDists(lat, longlat = longlat)
    rownames(lat_dists) <- colnames(lat_dists) <- names

    lon <- cbind(grid[, 1], mean(grid[, 2]))
    lon_dists <- sp::spDists(lon, longlat = longlat)
    rownames(lon_dists) <- colnames(lon_dists) <- names

    h <- sqrt(lon_dists^2 + lat_dists^2)
    rownames(h) <- colnames(h) <- names

    h1 <- matrix(0, ncol = n_var, nrow = n_var)
    for (i in 1:n_var) {
        for (j in 1:n_var) {
            h1[i, j] <- sign(grid[, 1][i] - grid[, 1][j]) *
                lon_dists[i, j]
        }
    }
    rownames(h1) <- colnames(h1) <- names

    h2 <- matrix(0, ncol = n_var, nrow = n_var)
    for (i in 1:n_var) {
        for (j in 1:n_var) {
            h2[i, j] <- sign(grid[, 2][i] - grid[, 2][j]) *
                lat_dists[i, j]
        }
    }
    rownames(h2) <- colnames(h2) <- names

    return(list(h = h, h1 = h1, h2 = h2))
}

#' Calculate (signed) distances between coordinates
#'
#' @inherit .find_dists return
#'
#' @param locations A matrix or data.frame of 2D points, the first column is
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
find_dists <- function(locations, longlat = TRUE) {
    if (NCOL(locations) != 2) {
        stop("`locations` must contain 2 columns", call. = FALSE)
    }

    names <- rownames(locations)
    if (any(table(names) > 1)) {
        stop("duplicate row names found in `locations`", call. = FALSE)
    }

    dists_ls <- .find_dists(locations, names = names, longlat = longlat)
    return(dists_ls)
}
