# Add lagr model outputted from [`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md) to a `mcgf` object.

Add lagr model outputted from
[`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md)
to a `mcgf` object.

## Usage

``` r
# S3 method for class 'mcgf'
add_lagr(x, fit_lagr, ...)

lagr(x) <- value
```

## Arguments

- x:

  An `mcgf` object.

- fit_lagr:

  Output from the
  [`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md)
  function.

- ...:

  Additional arguments. Not in use.

- value:

  A list containing the lagr model as well as its parameters. It must
  contains `model`, `par_lagr`, and `cor_lagr`.

## Value

`x` with newly added attributes of the Lagrangian model.

## See also

Other functions on fitting an mcgf:
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md),
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md),
[`fit_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf.md),
[`krige.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/krige.mcgf.md),
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
model(sim1_mcgf, old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Base-old model: 
#> - Parameters:
#> NULL
#> 
#> - Fixed parameters:
#> NULL
#> 
#> - Parameter estimation method: 
#> 
#> - Optimization function: 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma           a       alpha      nugget 
#> 0.001139264 0.500000000 0.627518829 0.734384140 0.000000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> ----------------------------------------
#>               Lagrangian
#> ----------------------------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>      lambda          v1          v2           k 
#>   0.1625903 224.5776224 188.5158088   2.0000000 
#> 
#> - Fixed parameters:
#> k 
#> 2 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
```
