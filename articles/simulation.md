# Simulating MCGF and RS-MCGF Processes

## Overview

`mcgf` can simulate Markov chain Gaussian fields (MCGF) and
regime-switching Markov chain Gaussian fields (RS-MCGF) directly from
user-specified covariance parameters. Simulation is useful for method
validation, sensitivity analysis, and learning the package workflow
before working with observed data.

This vignette focuses on the arguments that are most important in
practice.

## Distances

Start by defining the spatial locations and their signed distance
matrices.

``` r

library(mcgf)

locations <- matrix(
    c(0, 0, 2, 0, 1, 1.5),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(paste0("S", 1:3), c("x", "y"))
)

d <- find_dists(locations, longlat = FALSE)
d
#> $h
#>          S1       S2       S3
#> S1 0.000000 2.000000 1.802776
#> S2 2.000000 0.000000 1.802776
#> S3 1.802776 1.802776 0.000000
#> 
#> $h1
#>    S1 S2 S3
#> S1  0 -2 -1
#> S2  2  0  1
#> S3  1 -1  0
#> 
#> $h2
#>     S1  S2   S3
#> S1 0.0 0.0 -1.5
#> S2 0.0 0.0 -1.5
#> S3 1.5 1.5  0.0
```

`d$h` is Euclidean distance, while `d$h1` and `d$h2` retain horizontal
and vertical signs. The signed distances are needed by the Lagrangian
models.

## Simulate an MCGF

We first define a separable base covariance model.

``` r

par_base <- list(
    par_s = list(
        nugget = 0.05,
        c = 0.1,
        gamma = 0.5
    ),
    par_t = list(
        a = 0.2,
        alpha = 0.5
    )
)
```

A base-only first-order MCGF can be simulated by setting
`lagrangian = "none"` and `lambda = 0`.

``` r

set.seed(1)

x_base <- mcgf_sim(
    N = 200,
    base = "sep",
    lagrangian = "none",
    par_base = par_base,
    par_lagr = NULL,
    lambda = 0,
    dists = d,
    lag = 1
)

dim(x_base)
#> [1] 202   3
```

The first `lag + horizon` rows are initialization rows. If you want a
sample containing only generated observations, remove them before
analysis.

``` r

lag <- 1
horizon <- 1

x_base_use <- x_base[-seq_len(lag + horizon), , drop = FALSE]
```

## Add a Lagrangian component

The Lagrangian term models directional space-time asymmetry. For
example:

``` r

par_lagr <- list(v1 = 2, v2 = 1, k = 3)
```

``` r

set.seed(2)

x_lagr <- mcgf_sim(
    N = 200,
    base = "sep",
    lagrangian = "lagr_tri",
    par_base = par_base,
    par_lagr = par_lagr,
    lambda = 0.2,
    dists = d,
    lag = 2
)

dim(x_lagr)
#> [1] 203   3
```

`lambda` controls the contribution of the Lagrangian term. A value close
to zero gives a model dominated by the symmetric base covariance; a
larger value places more weight on directional advection.

## Simulate an RS-MCGF

For a regime-switching process, first supply a regime label for every
generated time point.

``` r

N <- 200
label <- rep(c(1, 2), each = N / 2)
table(label)
#> label
#>   1   2 
#> 100 100
```

Model arguments ending in `_ls` are lists. A length-one list can often
be reused across regimes; use one element per regime when parameters
differ.

``` r

par_base_1 <- par_base
par_base_2 <- list(
    par_s = list(
        nugget = 0.05,
        c = 0.2,
        gamma = 0.5
    ),
    par_t = list(
        a = 0.5,
        alpha = 0.7
    )
)
```

``` r

set.seed(3)

x_rs <- mcgf_rs_sim(
    N = N,
    label = label,
    base_ls = list("sep"),
    lagrangian_ls = list("none"),
    par_base_ls = list(par_base_1, par_base_2),
    par_lagr_ls = list(NULL),
    lambda_ls = list(0),
    dists_ls = list(d),
    lag_ls = list(1, 1)
)

head(x_rs)
#>   regime          S1          S2         S3
#> 1     NA  0.00000000  0.00000000  0.0000000
#> 2     NA  0.00000000  0.00000000  0.0000000
#> 3      1 -0.53172870 -0.51520683 -0.3898575
#> 4      1 -1.07997211 -0.85666860 -0.7907469
#> 5      1 -0.85276033 -0.28922534 -0.8124918
#> 6      1 -0.01006919  0.04511663 -0.5747100
```

The first column contains the regime label and the remaining columns
contain the simulated spatial process.

## Returning covariance details

Set `return_all = TRUE` when you also want the covariance matrices and
lag arrays used internally:

``` r

set.seed(4)

out <- mcgf_sim(
    N = 50,
    base = "sep",
    lagrangian = "none",
    par_base = par_base,
    par_lagr = NULL,
    lambda = 0,
    dists = d,
    lag = 1,
    return_all = TRUE
)

names(out)
#> [1] "X"   "par"
names(out$par)
#> [1] "cov_mat" "dists"   "u"
```

This is particularly useful when checking a simulation design or
comparing the implied covariance with an empirical estimate.

## Common choices

- `lag` determines the Markov order.
- `horizon` controls how many steps are generated jointly at each
  update.
- `sd` controls marginal standard deviations by location.
- `scale_time` rescales temporal lags before evaluating the covariance.
- `mu_c` and `mu_p` allow nonzero conditional means.
- `init` can supply explicit starting values.

For the covariance functions themselves, see
[`vignette("correlation-models", package = "mcgf")`](https://tianxia-jia.github.io/mcgf/articles/correlation-models.md).

## References
