# Add base model outputted from [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md) to an `mcgf_rs` object.

Add base model outputted from
[`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
to an `mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
add_base(x, fit_base_ls, fit_s_ls, fit_t_ls, sep = FALSE, old = FALSE, ...)
```

## Arguments

- x:

  An `mcgf_rs` object.

- fit_base_ls:

  Output from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function.

- fit_s_ls:

  Pure spatial model outputted from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function. Used only when `sep = TRUE`.

- fit_t_ls:

  Pure temporal model outputted from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function. Used only when `sep = TRUE`.

- sep:

  Logical; TRUE if spatial and temporal models are fitted separately.

- old:

  Logical; TRUE if the old base model needs to be kept. The lag and
  horizon of the new model are assumed to be the same as that of the old
  model.

- ...:

  Additional arguments. Not in use.

## Value

`x` with newly added attributes of the base model.

## Details

This function is equivalent to
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md)
for `mcgf_rs` objects.

After fitting the base model by
[`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md),
the results can be added to `x` by
[`add_base()`](https://tianxia-jia.github.io/mcgf/reference/add_base.md).
To supply the base model directly, use `base<-` to add the base model;
the value must contain the same output as
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md)
or `add_base.mcgf_rs()`.

## See also

Other functions on fitting an mcgf_rs:
[`add_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf_rs.md),
[`fit_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf_rs.md),
[`fit_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf_rs.md),
[`krige.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige.mcgf_rs.md),
[`krige_new.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/krige_new.mcgf_rs.md)

## Examples

``` r
data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
sim2_mcgf <- add_ccfs(sim2_mcgf, lag_max = 5)

# Fit a regime-switching pure spatial model
fit_spatial <-
    fit_base(
        sim2_mcgf,
        lag_ls = 5,
        model_ls = "spatial",
        par_init_ls = list(c(c = 0.00005, gamma = 0.5)),
        par_fixed_ls = list(c(nugget = 0))
    )

# Fit a regime-switching pure temporal model
fit_temporal <-
    fit_base(
        sim2_mcgf,
        lag_ls = 5,
        model_ls = "temporal",
        par_init_ls = list(
            list(a = 0.8, alpha = 0.8),
            list(a = 0.1, alpha = 0.1)
        )
    )

# Store the fitted models to 'sim2_mcgf'
sim2_mcgf <- add_base(sim2_mcgf,
    fit_s_ls = fit_spatial,
    fit_t_ls = fit_temporal,
    sep = TRUE
)

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

# Store the newly fitted model, and keep the old fit
sim2_mcgf <- add_base(sim2_mcgf, fit_base_ls = fit_sep, old = TRUE)
model(sim2_mcgf, model = "base", old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5, 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Regime switching: spatial: TRUE, temporal: TRUE 
#> --------------------
#>       Regime 1 
#> --------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma      nugget           a       alpha 
#> 0.001496696 0.500000000 0.000000000 0.437502762 0.302029968 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> --------------------
#>       Regime 2 
#> --------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma      nugget           a       alpha 
#> 0.004074408 0.500000000 0.000000000 0.330906238 0.874508851 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Regime switching: TRUE 
#> --------------------
#>       Regime 1 
#> --------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma           a       alpha      nugget 
#> 0.001539941 0.500000000 0.452841528 0.345950120 0.000000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> --------------------
#>       Regime 2 
#> --------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma           a       alpha      nugget 
#> 0.004107168 0.500000000 0.328388489 0.853896077 0.000000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
```
