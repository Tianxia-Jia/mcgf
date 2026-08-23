# Obtain kriging forecasts for new locations for an `mcgf_rs` object.

Obtain kriging forecasts for new locations for an `mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
krige_new(
  x,
  newdata = NULL,
  locations_new = NULL,
  dists_new_ls = NULL,
  newdata_new = NULL,
  sds_new_ls = 1,
  newlabel,
  soft = FALSE,
  prob,
  dists_new_base,
  model = c("all", "base"),
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

- locations_new:

  A matrix of data.frame of 2D points of new locations, first column
  longitude, second column latitude, both in decimal degrees. Supply
  only if `x` contains `locations`. Required when `dists_new_ls` is not
  supplied.

- dists_new_ls:

  List of signed distance matrices (vectors) with names `h`, `h1`, and
  'h2' for all locations(and for each regime), with new locations in the
  end. Each matrix must have the same number of columns. Required when
  `locations_new` is not supplied.

- newdata_new:

  Optional; a data.frame with the same number of rows as `newdata`. It
  contains the data of the new locations.

- sds_new_ls:

  List of the standard deviations of the new locations for each regime.
  Format must be the same as the output from
  [`sds.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/sds.mcgf.md).
  Default is 1 for all regimes.

- newlabel:

  A vector of new regime labels.

- soft:

  Logical; if true, soft forecasts (and bounds) are produced.

- prob:

  Matrix with simplex rows. Number of columns must be the same as unique
  labels in `x`.

- dists_new_base:

  Optional, list of distance matrices for the base model. Used when the
  base model is non-regime switching. Default is `h` from the first list
  of `dists_new_ls`.

- model:

  Which model to use. One of `all`, `base`, or `empirical`.

- interval:

  Logical; if TRUE, prediction intervals are computed.

- level:

  A numeric scalar between 0 and 1 giving the confidence level for the
  intervals (if any) to be calculated. Used when `interval = TRUE`

- ...:

  Additional arguments.

## Value

A list of kriging forecasts (and prediction intervals) for all
locations.

## Details

It produces simple kriging forecasts for a zero-mean mcgf for new
locations given their coordinates or relative distances. It supports
kriging for the `base` model and the `all` model which is the general
stationary model with the base and Lagrangian model from `x`.

Users can either supply the coordinates via `locations_new`, or a list
of distances for all locations via `dists_new_ls`, with new locations at
the end. `dists_new_ls` will be used to calculate the new covariance
matrices. When `locations_new` is used, make sure `x` contains the
attribute `locations` of the coordinates of the old locations. When
`dists_new_ls` is used, it should be a list of a list of signed distance
matrices of the same dimension, where each row corresponds to the
relative distances between a new location and old locations in the same
order as they appear in `x`. If only one list is provided, it will be
used for all regimes.

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
[`krige.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige.mcgf_rs.md)

## Examples

``` r
data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data,
    locations = sim2$locations,
    label = sim2$label
)
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

# Calculate the simple kriging predictions and intervals for all locations
locations_new <- rbind(c(-110, 55), c(-109, 54))
sim2_krige <- krige_new(sim2_mcgf,
    locations_new = locations_new,
    model = "base", interval = TRUE
)
#> Latitude for new locations is outside of the grid.
```
