# Obtain kriging forecasts for new locations for an `mcgf` object.

Obtain kriging forecasts for new locations for an `mcgf` object.

## Usage

``` r
# S3 method for class 'mcgf'
krige_new(
  x,
  newdata = NULL,
  locations_new = NULL,
  dists_new = NULL,
  newdata_new = NULL,
  sds_new = 1,
  model = c("all", "base"),
  interval = FALSE,
  level = 0.95,
  dists_new_base,
  ...
)
```

## Arguments

- x:

  An `mcgf` object.

- newdata:

  A data.frame with the same column names as `x`. If `newdata` is
  missing the forecasts at the original data points are returned.

- locations_new:

  A matrix of data.frame of 2D points of new locations, first column
  longitude, second column latitude, both in decimal degrees. Supply
  only if `x` contains `locations`. Required when `dists_new` is not
  supplied.

- dists_new:

  List of signed distance matrices (vectors) with names `h`, `h1`, and
  'h2' for all locations, with new locations in the end. Each matrix
  must have the same number of columns. Required when `locations_new` is
  not supplied.

- newdata_new:

  Optional; a data.frame with the same number of rows as `newdata`. It
  contains the data of the new locations.

- sds_new:

  The standard deviations of the new locations. Default is 1.

- model:

  Which model to use. One of `all` or `base`.

- interval:

  Logical; if TRUE, prediction intervals are computed.

- level:

  A numeric scalar between 0 and 1 giving the confidence level for the
  intervals (if any) to be calculated. Used when `interval = TRUE`.

- dists_new_base:

  Optional; a distance array for all locations for the base model, with
  new locations in the end. Used for the base model.

- ...:

  Additional arguments. Not in use.

## Value

A list of kriging forecasts (and prediction intervals) for all
locations.

## Details

It produces simple kriging forecasts for a zero-mean mcgf for new
locations given their coordinates or relative distances. It supports
kriging for the `base` model and the `all` model which is the general
stationary model with the base and Lagrangian model from `x`.

Users can either supply the coordinates via `locations_new`, or a list
of distances for all locations via `dists_new`, with new locations at
the end. `dists_new` will be used to calculate the new covariance
matrices. When `locations_new` is used, make sure `x` contains the
attribute `locations` of the coordinates of the old locations. When
`dists_new` is used, it should be a list of signed distance matrices of
the same dimension, where each row corresponds to the relative distances
between a new location and old locations in the same order as they
appear in `x`.

If data for the new locations are available, use `newdata_new` to
include them and they will be used to calculate the kriging forecasts
for the new locations; otherwise only data of the old locations will be
used via `newdata`.

When `interval = TRUE`, a confidence interval for each forecasts and
each horizon is given. Note that it does not compute confidence regions.

## See also

Other functions on fitting an mcgf:
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md),
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md),
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md),
[`fit_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf.md),
[`krige.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/krige.mcgf.md)

## Examples

``` r
data(sim1)
sim1_mcgf <- mcgf(sim1$data, locations = sim1$locations)
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

# Calculate the simple kriging predictions and intervals for all locations
locations_new <- rbind(c(-110, 55), c(-109, 54))
sim1_krige <- krige_new(sim1_mcgf,
    locations_new = locations_new,
    interval = TRUE
)
#> Latitude for new locations is outside of the grid.
```
