new_mcgf <- function(data, locations, dists) {

    data <- as.data.frame(data)

    if (!missing(locations)) {
        dists <- find_dists(locations)
    }
        structure(.Data = data, dists = dists, class = c("data.frame", "mcgf"))
}

validate_mcgf <- function(x) {
    x
}

mcgf <- function(data, locations, dists, time) {

    if (!is.data.frame(data) && !is.matrix(data))
        stop("'data' must be a matrix or data.frame.", call. = FALSE)

    if (missing(locations) && missing(dists)) {
        stop("either 'locations' or 'dists' needs to be provided.",
             call. = FALSE)
    }

    if (!missing(locations) && !missing(dists)) {
        stop("do not provide both 'locations' and 'dists'.", call. = FALSE)
    }

    name_var <- colnames(data)

    n_var <- NCOL(data)

    if (missing(time))
        warning("'time' not provided, row names of 'data' are used.")
        time <- 1:NROW(data)

    if (length(time) != NROW(data)) {
        stop("length of 'time' must be the same as the number of rows of ",
             "'data'.", call. = FALSE)
    }

    diff_time <- diff(time)
    if (length(unique(diff_time)) != 1)
        stop("time must be equally spaced.")
    if (unique(diff_time) < 0)
        stop("time must be in ascending order.")

    rownames(data) <- time

    if (!missing(locations)) {
        if (!missing(dists)) {
            stop("do not provide both 'locations' and 'dists'.", call. = FALSE)
        } else {
            if (ncol(data) != nrow(locations))
                stop("number of columns of 'data' must be the same as the ",
                     "number of rows of 'stations'", call. = FALSE)
            if (any(colnames(data) != rownames(locations))) {
                warning("rownames of 'locations' not the same as the column ",
                        "names of 'data', the latter is used.", call. = FALSE)
                rownames(locations) <- colnames(data)
            }
            return(validate_mcgf(new_mcgf(data = data, locations = locations)))
        }
    } else {
        if (!is.list(dists)) {
            stop("'dists' must be a list.", call. = FALSE)
        }
        if (any(!c("h1", "h2") %in% names(dists)))
            stop("'dists' must contain 'h1' and 'h2',", call. = FALSE)

        if (!is.matrix(dists$h1))
            stop("'h1' in 'dists' must be a matrix.", call. = FALSE)
        if (!is.matrix(dists$h2))
            stop("'h2' in 'dists' must be a matrix.", call. = FALSE)

        if (any(dim(dists$h1) != c(n_var, n_var)))
            stop("'h1' in 'dists' must be a matrix of dimension ",
                 n_var, " x ", n_var, ".", call. = FALSE)
        if (any(dim(dists$h2) != c(n_var, n_var)))
            stop("'h2' in 'dists' must be a matrix of dimension ",
                 n_var, " x ", n_var, ".", call. = FALSE)

        check_dist_sign(dists$h1, name = "h1")
        check_dist_sign(dists$h2, name = "h2")

        if (is.null(rownames(dists$h1)))
            rownames(dists$h1) <- colnames(data)
        if (is.null(colnames(dists$h1)))
            colnames(dists$h1) <- colnames(data)
        if (is.null(rownames(dists$h2)))
            rownames(dists$h2) <- colnames(data)
        if (is.null(colnames(dists$h2)))
            colnames(dists$h2) <- colnames(data)

        if (!is.null(dists$h)) {
            if (!is.matrix(dists$h))
                stop("'h' in 'dists' must be a matrix.", call. = FALSE)

            if (any(dim(dists$h) != c(n_var, n_var)))
                stop("'h' in 'dists' must be a matrix of dimension ",
                     n_var, " x ", n_var, ".", call. = FALSE)

            check_dist(x = dists$h, "h")

            if (is.null(rownames(dists$h)))
                rownames(dists$h) <- colnames(data)
            if (is.null(colnames(dists$h)))
                colnames(dists$h) <- colnames(data)

        } else {
            dists$h <- sqrt(dists$h1 ^ 2 + dists$h2 ^ 2)
        }
        return(validate_mcgf(new_mcgf(data = data, dists = dists)))
    }
}
