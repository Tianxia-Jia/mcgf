cor2cov <- function(cov, sd) {

    stopifnot(dim(cov) == c(length(sd), length(sd)))

    sd_mat <- as.matrix(sd)
    cov * sd_mat %*% t(sd_mat)
}
