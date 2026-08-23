# Parameter estimation for Lagrangian correlation functions for an `mcgf` object.

Parameter estimation for Lagrangian correlation functions for an `mcgf`
object.

## Usage

``` r
# S3 method for class 'mcgf'
fit_lagr(
  x,
  model = c("lagr_tri", "lagr_askey", "lagr_exp", "none"),
  method = c("wls", "mle"),
  optim_fn = c("nlminb", "optim", "other"),
  par_fixed = NULL,
  par_init,
  lower = NULL,
  upper = NULL,
  other_optim_fn = NULL,
  dists_base = FALSE,
  dists_lagr = NULL,
  ...
)
```

## Arguments

- x:

  An `mcgf` object containing attributes `dists`, `acfs`, `ccfs`, and
  `sds`. `x` must have been passed to
  [`add_base()`](https://tianxia-jia.github.io/mcgf/reference/add_base.md)
  or `base<-`

- model:

  Base model, one of `lagr_tri`, `lagr_askey`, `lagr_exp`, or `none`. If
  `none`, NULLs are returned.

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

  Optional; lower bounds of parameters lambda, v1, v2, and k.

- upper:

  Optional: upper bounds of parameters lambda, v1, v2, and k.

- other_optim_fn:

  Optional, other optimization functions. The first two arguments must
  be initial values for the parameters and a function to be minimized
  respectively (same as that of `optim` and `nlminb`).

- dists_base:

  Logical; if TRUE `dists_base` from the base model is used as the
  distance.

- dists_lagr:

  List of distance matrices/arrays. Used when `dists_base` is FALSE. If
  NULL, `dists(x)` is used.

- ...:

  Additional arguments passed to `optim_fn`.

## Value

A list containing outputs from optimization functions of `optim_fn`.

## Details

This function fits the Lagrangian models using weighted least squares
and maximum likelihood estimation. The base model must be fitted first
using
[`add_base()`](https://tianxia-jia.github.io/mcgf/reference/add_base.md)
or `base<-`. Optimization functions such as `nlminb` and `optim` are
supported. Other functions are also supported by setting
`optim_fn = "other"` and supplying `other_optim_fn`. `lower` and `upper`
are lower and upper bounds of parameters in `par_init` and default
bounds are used if they are not specified.

Note that both `wls` and `mle` are heuristic approaches when `x`
contains observations from a subset of the discrete spatial domain,
though estimation results are close to that using the full spatial
domain for large sample sizes.

Since parameters for the base model and the Lagrangian model are
estimated sequentially, more accurate estimation may be obtained if the
full model is fitted all at once.

## See also

Other functions on fitting an mcgf:
[`add_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf.md),
[`add_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf.md),
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md),
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
fit_lagr$fit
#> $par
#>      lambda          v1          v2 
#>   0.1625903 224.5776224 188.5158088 
#> 
#> $objective
#> [1] 1.533448
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 34
#> 
#> $evaluations
#> function gradient 
#>       37      118 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
```
