test_that("cor_lagr_tri() matches its defining formula", {
    h1 <- matrix(c(0, -1, 1, 0), 2)
    h2 <- matrix(c(0, -2, 2, 0), 2)
    u <- matrix(0.25, 2, 2)
    v1 <- 2
    v2 <- 1
    k <- 3

    v_norm <- sqrt(v1^2 + v2^2)
    expected <- 1 - abs((h1 * v1 + h2 * v2) / v_norm - v_norm * u) /
        (k * v_norm)
    expected[expected < 0] <- 0

    expect_equal(cor_lagr_tri(v1, v2, k, h1, h2, u), expected)
})

test_that("cor_lagr_askey() matches its defining formula", {
    h1 <- matrix(c(0, -1, 1, 0), 2)
    h2 <- matrix(c(0, -2, 2, 0), 2)
    u <- matrix(0.25, 2, 2)
    v1 <- 2
    v2 <- 1
    k <- 3

    v_norm <- sqrt(v1^2 + v2^2)
    r <- 1 - sqrt((h1 - v1 * u)^2 + (h2 - v2 * u)^2) /
        (k * v_norm)
    r[r < 0] <- 0

    expect_equal(cor_lagr_askey(v1, v2, k, h1, h2, u), r^(3 / 2))
})

test_that("cor_lagr_exp() matches its defining formula", {
    h1 <- matrix(c(0, -1, 1, 0), 2)
    h2 <- matrix(c(0, -2, 2, 0), 2)
    u <- matrix(0.25, 2, 2)
    v1 <- 2
    v2 <- 1
    k <- 3

    v_norm <- sqrt(v1^2 + v2^2)
    expected <- exp(-sqrt((h1 - v1 * u)^2 + (h2 - v2 * u)^2) / (k * v_norm))

    expect_equal(cor_lagr_exp(v1, v2, k, h1, h2, u), expected)
})

test_that("Lagrangian functions validate dimensions and parameters", {
    h1 <- matrix(c(0, -1, 1, 0), 2)
    h2 <- matrix(c(0, -2, 2, 0), 2)
    u <- matrix(0, 2, 2)

    for (fn in list(cor_lagr_tri, cor_lagr_askey, cor_lagr_exp)) {
        expect_error(fn(0, 0, 2, h1, h2, u), "cannot be zero")
        expect_error(fn(1, 1, 0, h1, h2, u), "`k` must be positive")
        expect_error(fn(1, 1, 2, h1, h2[1, , drop = FALSE], u))
        expect_error(fn(1, 1, 2, h1, h2, array(c(2, 2, 2))))
    }
})
