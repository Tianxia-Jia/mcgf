# Forecasting at New Locations

## Overview

[`krige_new()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.md)
extends an estimated MCGF model to spatial locations that were not used
for fitting. This is useful when:

- a station is newly installed;
- forecasts are required on a target grid;
- some sites are intentionally held out for validation; or
- future scenario generation is needed at additional locations.

There are two ways to describe the new sites:

1.  provide their coordinates with `locations_new`, or
2.  provide precomputed signed distances with `dists_new`.

This vignette demonstrates the coordinate workflow.

## Prepare an MCGF object

The built-in `sim1` data contain both observations and locations.

``` r

library(mcgf)
data(sim1)

x <- mcgf(
    sim1$data,
    locations = sim1$locations
)
#> `time` is not provided, assuming rows are equally spaced temporally.

x <- add_acfs(x, lag_max = 5)
x <- add_ccfs(x, lag_max = 5)
```

## Fit the base covariance model

For a first new-location forecast, a separable base model is sufficient.

``` r

fit_sep <- fit_base(
    x,
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

x <- add_base(x, fit_base = fit_sep)
```

## Define a new location

For illustration, place a new site near the centre of the existing
network.

``` r

new_location <- matrix(
    colMeans(sim1$locations),
    nrow = 1
)
rownames(new_location) <- "New_site"

new_location
#>               [,1]     [,2]
#> New_site -109.9254 50.20862
```

## Forecast with `krige_new()`

``` r

pred <- krige_new(
    x,
    locations_new = new_location,
    model = "base",
    interval = TRUE
)

dim(pred$fit)
#> [1] 1000   11
```

The output contains forecasts for the original sites **and** the
appended new site. The new site is the last spatial column.

``` r

new_index <- dim(pred$fit)[2]

head(pred$fit[, new_index])
#>          1          2          3          4          5          6 
#>         NA         NA         NA         NA         NA 0.01776882
head(pred$lower[, new_index])
#>         1         2         3         4         5         6 
#>        NA        NA        NA        NA        NA -1.534888
head(pred$upper[, new_index])
#>        1        2        3        4        5        6 
#>       NA       NA       NA       NA       NA 1.570425
```

## Forecast using new observed data

If you have a separate future period for the original locations, supply
it with `newdata`:

``` r

pred_future <- krige_new(
    x,
    newdata = future_data,
    locations_new = new_location,
    model = "base",
    interval = TRUE
)
```

`future_data` must have the same columns, in the same order, as the
fitted `mcgf` object.

## If observations at the new location are available

You can also provide `newdata_new`. These observations are then included
in the lagged information used for forecasting.

``` r

pred_with_new_history <- krige_new(
    x,
    newdata = future_data,
    locations_new = new_location,
    newdata_new = future_new_site,
    model = "base",
    interval = TRUE
)
```

The number of rows in `newdata_new` must match `newdata`, and the number
of columns must match the number of new locations.

## Supplying distances directly

If the object does not contain coordinates, calculate the distances
yourself and pass them through `dists_new`.

``` r

d_new <- find_dists_new(
    locations = observed_locations,
    locations_new = new_location,
    longlat = TRUE
)

pred <- krige_new(
    x,
    dists_new = d_new,
    model = "base"
)
```

Do not supply both `locations_new` and `dists_new`.

## General stationary model

Once a Lagrangian model has been fitted and added with
[`add_lagr()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.md),
switch to:

``` r

pred_all <- krige_new(
    x,
    locations_new = new_location,
    model = "all",
    interval = TRUE
)
```

This extends both the base and directional Lagrangian covariance to the
new location.

## Regime-switching objects

[`krige_new()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.md)
also has an `mcgf_rs` method. Use `dists_new_ls` and `sds_new_ls` when
distance or marginal-variance assumptions differ by regime. For soft
regime forecasts, supply regime probabilities using `prob`.

See
[`vignette("mcgf_rs", package = "mcgf")`](https://tianxia-jia.github.io/mcgf/articles/mcgf_rs.md)
for a complete RS-MCGF fit.
