# Parameter estimation for symmetric correlation functions for an `mcgf` object.

Parameter estimation for symmetric correlation functions for an `mcgf`
object.

## Usage

``` r
# S3 method for class 'mcgf'
fit_base(
  x,
  lag,
  horizon = 1,
  model = c("spatial", "temporal", "sep", "fs", "none"),
  method = c("wls", "mle"),
  optim_fn = c("nlminb", "optim", "other"),
  par_fixed = NULL,
  par_init,
  lower = NULL,
  upper = NULL,
  other_optim_fn = NULL,
  dists_base = NULL,
  scale_time = 1,
  ...
)
```

## Arguments

- x:

  An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
  `sds`.

- lag:

  Integer time lag.

- horizon:

  Integer forecast horizon.

- model:

  Base model, one of `spatial`, `temporal`, `sep`, `fs`, `none`. Only
  `sep` and `fs` are supported when `method = mle`. If `none`, NULLs are
  returned.

- method:

  Parameter estimation methods, weighted least square (`wls`) or maximum
  likelihood estimation (`mle`).

- optim_fn:

  Optimization functions, one of `nlminb`, `optim`, `other`. When
  `optim_fn = other`, supply `other_optim_fn`.

- par_fixed:

  Fixed parameters.

- par_init:

  Initial values for parameters to be optimized.

- lower:

  Optional; lower bounds of parameters.

- upper:

  Optional: upper bounds of parameters.

- other_optim_fn:

  Optional, other optimization functions. The first two arguments must
  be initial values for the parameters and a function to be minimized
  respectively (same as that of `optim` and `nlminb`).

- dists_base:

  List of distance matrices. If NULL, `dists(x)` is used. Must be a
  matrix or an array of distance matrices.

- scale_time:

  Scale of time unit, default is 1. `lag` is divided by `scale_time` for
  parameter estimation.

- ...:

  Additional arguments passed to `optim_fn`.

## Value

A list containing outputs from optimization functions of `optim_fn`.

## Details

This function fits the separable and fully symmetric models using
weighted least squares and maximum likelihood estimation. Optimization
functions such as `nlminb` and `optim` are supported. Other functions
are also supported by setting `optim_fn = "other"` and supplying
`other_optim_fn`. `lower` and `upper` are lower and upper bounds of
parameters in `par_init` and default bounds are used if they are not
specified.

Note that both `wls` and `mle` are heuristic approaches when `x`
contains observations from a subset of the discrete spatial domain,
though estimation results are close to that using the full spatial
domain for large sample sizes.

## See also

Other functions on fitting an mcgf:
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md),
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md),
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
fit_spatial$fit
#> $par
#>           c       gamma 
#> 0.001145201 0.500000000 
#> 
#> $objective
#> [1] 1.385845
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 7
#> 
#> $evaluations
#> function gradient 
#>       18       17 
#> 
#> $message
#> [1] "both X-convergence and relative convergence (5)"
#> 

# Fit a pure temporal model
fit_temporal <- fit_base(
    sim1_mcgf,
    model = "temporal",
    lag = 5,
    par_init = c(a = 0.3, alpha = 0.5)
)
fit_temporal$fit
#> $par
#>         a     alpha 
#> 0.6534102 0.7572530 
#> 
#> $objective
#> [1] 0.004445817
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 16
#> 
#> $evaluations
#> function gradient 
#>       18       38 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 

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
fit_sep$fit
#> $par
#>           c       gamma           a       alpha 
#> 0.001139264 0.500000000 0.627518829 0.734384140 
#> 
#> $objective
#> [1] 3.196846
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 25
#> 
#> $evaluations
#> function gradient 
#>       46      113 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
```
