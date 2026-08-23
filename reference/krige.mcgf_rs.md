# Obtain kriging forecasts for an `mcgf_rs` object.

Obtain kriging forecasts for an `mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
krige(
  x,
  newdata = NULL,
  newlabel = NULL,
  soft = FALSE,
  prob,
  model = c("all", "base", "empirical"),
  interval = FALSE,
  level = 0.95,
  ...
)
```

## Arguments

- x:

  An `mcgf_rs` object.

- newdata:

  A data.frame with the same column names as `x`. If `newdata` is
  missing the forecasts at the original data points are returned.

- newlabel:

  A vector of new regime labels.

- soft:

  Logical; if true, soft forecasts (and bounds) are produced.

- prob:

  Matrix with simplex rows. Number of columns must be the same as unique
  labels in `x`.

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

When `soft = TRUE`, `prob` will be used to compute the soft forecasts
(weighted forecasts). The number of columns must match the number of
unique levels in `x`. The column order must be the same as the order of
regimes as in `levels(attr(x, "label", exact = TRUE))`. If not all
regimes are seen in `newlabel`, then only relevant columns in `prob` are
used.

When `interval = TRUE`, confidence interval for each forecasts and each
horizon is given. Note that it does not compute confidence regions.

## See also

Other functions on fitting an mcgf_rs:
[`add_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf_rs.md),
[`add_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf_rs.md),
[`fit_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf_rs.md),
[`fit_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf_rs.md),
[`krige_new.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.mcgf_rs.md)

## Examples

``` r
data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)

# Fit a regime-switching separable model
fit_sep <- fit_base(
    sim2_mcgf,
    lag_ls = 5,
    model_ls = "sep",
    par_init_ls = list(list(
        c = 0.00005,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.5
    )),
    par_fixed_ls = list(c(nugget = 0))
)

# Store the fitted separable models to 'sim2_mcgf'
sim2_mcgf <- add_base(sim2_mcgf, fit_base_ls = fit_sep)

# Calculate the simple kriging predictions and intervals
sim2_krige <- krige(sim2_mcgf, model = "base", interval = TRUE)

# Calculate RMSE for each location
rmse <- sqrt(colMeans((sim2_mcgf - sim2_krige$fit)^2, na.rm = TRUE))
rmse
#>         1         2         3         4         5         6         7         8 
#> 0.6200119 0.6352447 0.6427659 0.6505784 0.6567752 0.6380914 0.6288325 0.6552253 
#>         9        10 
#> 0.6387179 0.6415569 

# Calculate MAE for each location
mae <- colMeans(abs(sim2_mcgf - sim2_krige$fit), na.rm = TRUE)
mae
#>         1         2         3         4         5         6         7         8 
#> 0.4849978 0.4971572 0.5105414 0.5119415 0.5208628 0.5077890 0.4947488 0.5237569 
#>         9        10 
#> 0.4990462 0.5048435 

# Calculate POPI for each location
popi <- colMeans(
    sim2_mcgf < sim2_krige$lower | sim2_mcgf > sim2_krige$upper,
    na.rm = TRUE
)
popi
#>          1          2          3          4          5          6          7 
#> 0.05427136 0.05025126 0.04824121 0.04422111 0.04824121 0.05527638 0.03819095 
#>          8          9         10 
#> 0.04221106 0.06231156 0.05829146 
```
