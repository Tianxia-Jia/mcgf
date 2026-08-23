# Add base model outputted from [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md) to an `mcgf` object.

Add base model outputted from
[`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
to an `mcgf` object.

## Usage

``` r
# S3 method for class 'mcgf'
add_base(x, fit_base, fit_s, fit_t, sep = FALSE, old = FALSE, ...)

base(x) <- value
```

## Arguments

- x:

  An `mcgf` object.

- fit_base:

  Output from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function.

- fit_s:

  Pure spatial model outputted from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function. Used only when `sep = TRUE`.

- fit_t:

  Pure temporal model outputted from the
  [`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md)
  function. Used only when `sep = TRUE`.

- sep:

  Logical; TRUE if spatial and temporal models are fitted separately.

- old:

  Logical; TRUE if the old base model needs to be kept.

- ...:

  Additional arguments. Not in use.

- value:

  A list containing the base model as well as its parameters. It must
  contains the same output as `add_base.mcgf()` or
  [`add_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf_rs.md).

## Value

`x` with newly added attributes of the base model.

## Details

After fitting the base model by
[`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md),
the results can be added to `x` by
[`add_base()`](https://tianxia-jia.github.io/mcgf/reference/add_base.md).
To supply the base model directly, use `base<-` to add the base model;
the value must contain `model`, `par_base`, `cor_base`, `lag`, and
`horizon`.

## See also

Other functions on fitting an mcgf:
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md),
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

# Fit a pure spatial model
fit_spatial <- fit_base(
    sim1_mcgf,
    model = "spatial",
    lag = 5,
    par_init = c(c = 0.001, gamma = 0.5),
    par_fixed = c(nugget = 0)
)
# Fit a pure temporal model
fit_temporal <- fit_base(
    sim1_mcgf,
    model = "temporal",
    lag = 5,
    par_init = c(a = 0.3, alpha = 0.5)
)

# Store the fitted models to 'sim1_mcgf'
sim1_mcgf <-
    add_base(sim1_mcgf,
        fit_s = fit_spatial,
        fit_t = fit_temporal,
        sep = TRUE
    )

# Fit a separable model
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
# Store the newly fitted model, and keep the old fit
sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep, old = TRUE)
model(sim1_mcgf, model = "base", old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Base-old model: sep 
#> - Parameters:
#>           c       gamma      nugget           a       alpha 
#> 0.001145201 0.500000000 0.000000000 0.653410220 0.757253028 
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
```
