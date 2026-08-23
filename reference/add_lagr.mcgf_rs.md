# Add lagr model outputted from [`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md) to a `mcgf_rs` object.

Add lagr model outputted from
[`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md)
to a `mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
add_lagr(x, fit_lagr_ls, ...)
```

## Arguments

- x:

  An `mcgf_rs` object.

- fit_lagr_ls:

  Output from the
  [`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md)
  function.

- ...:

  Additional arguments. Not in use.

## Value

`x` with newly added attributes of the Lagrangian model.

## Details

This function is equivalent to
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md)
for `mcgf_rs` objects.

After fitting the Lagrangian model by
[`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md),
the results can be added to `x` by
[`add_lagr()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.md).
To supply the Lagrangian model directly, use `lagr<-` to add the
Lagrangian model; the value must contain the same output as
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md)
or `add_lagr.mcgf_rs()`.

## See also

Other functions on fitting an mcgf_rs:
[`add_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf_rs.md),
[`fit_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf_rs.md),
[`fit_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf_rs.md),
[`krige.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige.mcgf_rs.md),
[`krige_new.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.mcgf_rs.md)

## Examples

``` r
data(sim3)
sim3_mcgf <- mcgf_rs(sim3$data, dists = sim3$dists, label = sim3$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim3_mcgf <- add_acfs(sim3_mcgf, lag_max = 5)
sim3_mcgf <- add_ccfs(sim3_mcgf, lag_max = 5)

# Fit a fully symmetric model with known variables
fit_fs <- fit_base(
    sim3_mcgf,
    lag_ls = 5,
    model_ls = "fs",
    rs = FALSE,
    par_init_ls = list(list(beta = 0)),
    par_fixed_ls = list(list(
        nugget = 0,
        c = 0.05,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.2
    ))
)

# Set beta to 0 to fit a separable model with known variables
fit_fs[[1]]$fit$par <- 0

# Store the fitted separable model to 'sim3_mcgf'
sim3_mcgf <- add_base(sim3_mcgf, fit_base_ls = fit_fs)

# Fit a regime-switching Lagrangian model.
fit_lagr_rs <- fit_lagr(
    sim3_mcgf,
    model_ls = list("lagr_tri"),
    par_init_ls = list(
        list(v1 = -50, v2 = 50),
        list(v1 = 100, v2 = 100)
    ),
    par_fixed_ls = list(list(lambda = 0.2, k = 2))
)

# Store the fitted Lagrangian model to 'sim3_mcgf'
sim3_mcgf <- add_lagr(sim3_mcgf, fit_lagr_ls = fit_lagr_rs)
model(sim3_mcgf)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5, 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Regime switching: FALSE
#> - Base model: fs 
#> - Parameters:
#>   beta nugget      c  gamma      a  alpha 
#>   0.00   0.00   0.05   0.50   0.50   0.20 
#> 
#> - Fixed parameters:
#> nugget      c  gamma      a  alpha 
#>   0.00   0.05   0.50   0.50   0.20 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> ----------------------------------------
#>               Lagrangian
#> ----------------------------------------
#> - Regime switching: TRUE 
#> --------------------
#>       Regime 1 
#> --------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>        v1        v2    lambda         k 
#> -106.5181  119.4434    0.2000    2.0000 
#> 
#> - Fixed parameters:
#> lambda      k 
#>    0.2    2.0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> --------------------
#>       Regime 2 
#> --------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>       v1       v2   lambda        k 
#> 210.1156 242.3082   0.2000   2.0000 
#> 
#> - Fixed parameters:
#> lambda      k 
#>    0.2    2.0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
```
