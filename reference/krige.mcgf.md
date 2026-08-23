# Obtain kriging forecasts for an `mcgf` object.

Obtain kriging forecasts for an `mcgf` object.

## Usage

``` r
# S3 method for class 'mcgf'
krige(
  x,
  newdata = NULL,
  model = c("all", "base", "empirical"),
  interval = FALSE,
  level = 0.95,
  ...
)
```

## Arguments

- x:

  An `mcgf` object.

- newdata:

  A data.frame with the same column names as `x`. If `newdata` is
  missing the forecasts at the original data points are returned.

- model:

  Which model to use. One of `all`, `base`, or `empirical`.

- interval:

  Logical; if TRUE, prediction intervals are computed.

- level:

  A numeric scalar between 0 and 1 giving the confidence level for the
  intervals (if any) to be calculated. Used when `interval = TRUE`

- ...:

  Additional arguments. Give `lag` and `horizon` if they are not defined
  in `x` for the `empirical` model.

## Value

A list of kriging forecasts (and prediction intervals).

## Details

It produces simple kriging forecasts for a zero-mean mcgf. It supports
kriging for the `empirical` model, the `base` model, and the `all` model
which is the general stationary model with the base and Lagrangian model
from `x`.

When `interval = TRUE`, a confidence interval for each forecasts and
each horizon is given. Note that it does not compute confidence regions.

## See also

Other functions on fitting an mcgf:
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md),
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md),
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md),
[`fit_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf.md),
[`krige_new.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.mcgf.md)

## Examples

``` r
data(sim1)
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = 5)

# Fit a separable model and store it to 'sim1_mcgf'
fit_sep <- fit_base(
    sim1_mcgf,
    model = "sep",
    lag = 5,
    par_init = c(
        c = 0.001,
        gamma = 0.5,
        a = 0.3,
        alpha = 0.5
    ),
    par_fixed = c(nugget = 0)
)
sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep)

# Fit a Lagrangian model
fit_lagr <- fit_lagr(
    sim1_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 300, v2 = 300, lambda = 0.15),
    par_fixed = c(k = 2)
)

# Store the fitted Lagrangian model to 'sim1_mcgf'
sim1_mcgf <- add_lagr(sim1_mcgf, fit_lagr = fit_lagr)

# Calculate the simple kriging predictions and intervals
sim1_krige <- krige(sim1_mcgf, interval = TRUE)

# Calculate RMSE for each location
rmse <- sqrt(colMeans((sim1_mcgf - sim1_krige$fit)^2, na.rm = TRUE))
rmse
#>         1         2         3         4         5         6         7         8 
#> 0.7303482 0.7327718 0.7356240 0.7345241 0.7334888 0.7359809 0.7336117 0.7330753 
#>         9        10 
#> 0.7363262 0.7446223 

# Calculate MAE for each location
mae <- colMeans(abs(sim1_mcgf - sim1_krige$fit), na.rm = TRUE)
mae
#>         1         2         3         4         5         6         7         8 
#> 0.5835321 0.5836252 0.5845320 0.5853083 0.5798051 0.5922677 0.5855238 0.5824134 
#>         9        10 
#> 0.5843573 0.5825537 

# Calculate POPI for each location
popi <- colMeans(
    sim1_mcgf < sim1_krige$lower | sim1_mcgf > sim1_krige$upper,
    na.rm = TRUE
)
popi
#>          1          2          3          4          5          6          7 
#> 0.04623116 0.04120603 0.05025126 0.04623116 0.05226131 0.04020101 0.04924623 
#>          8          9         10 
#> 0.04321608 0.04020101 0.04422111 
```
