# Parameter estimation for symmetric correlation functions for an `mcgf_rs` object.

Parameter estimation for symmetric correlation functions for an
`mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
fit_base(
  x,
  lag_ls,
  horizon = 1,
  model_ls,
  method_ls = "wls",
  optim_fn_ls = "nlminb",
  par_fixed_ls = list(NULL),
  par_init_ls,
  lower_ls = list(NULL),
  upper_ls = list(NULL),
  other_optim_fn_ls = list(NULL),
  dists_base_ls = list(NULL),
  scale_time = 1,
  rs = TRUE,
  ...
)
```

## Arguments

- x:

  An `mcgf_rs` object containing attributes `dists`, `acfs`, `ccfs`, and
  `sds`.

- lag_ls:

  List of integer time lags.

- horizon:

  Integer forecast horizon.

- model_ls:

  List of base models, each element must be one of `spatial`,
  `temporal`, `sep`, `fs`. Only `sep` and `fs` are supported when `mle`
  is used in `model_ls`.

- method_ls:

  List of parameter estimation methods, weighted least square (`wls`) or
  maximum likelihood estimation (`mle`).

- optim_fn_ls:

  List of optimization functions, each element must be one of `nlminb`,
  `optim`, `other`. When use `other`, supply `other_optim_fn_ls`.

- par_fixed_ls:

  List of fixed parameters.

- par_init_ls:

  List of initial values for parameters to be optimized.

- lower_ls:

  Optional; list of lower bounds of parameters.

- upper_ls:

  Optional: list of upper bounds of parameters.

- other_optim_fn_ls:

  Optional, list of other optimization functions. The first two
  arguments must be initial values for the parameters and a function to
  be minimized respectively (same as that of `optim` and `nlminb`).

- dists_base_ls:

  List of lists of distance matrices. If NULL, `dists(x)` is used. Each
  element must be a matrix or an array of distance matrices.

- scale_time:

  Scale of time unit, default is 1. `lag` is divided by `scale_time` for
  parameter estimation.

- rs:

  Logical; if TRUE `x` is treated as a regime-switching, FALSE if the
  parameters need to be estimated in a non-regime-switching setting.

- ...:

  Additional arguments passed to all `optim_fn_ls`.

## Value

A list containing outputs from optimization functions of `optim_fn` for
each regime.

## Details

This functions is the regime-switching variant of
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md).
Arguments are in lists. The length of arguments that end in `_ls` must
be 1 or the same as the number of regimes in `x`. If the length of an
argument is 1, then it is set the same for all regimes. Refer to
[`fit_base.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf.md)
for more details of the arguments.

Note that both `wls` and `mle` are heuristic approaches when `x`
contains observations from a subset of the discrete spatial domain,
though estimation results are close to that using the full spatial
domain for large sample sizes.

## See also

Other functions on fitting an mcgf_rs:
[`add_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf_rs.md),
[`add_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf_rs.md),
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
lapply(fit_spatial[1:2], function(x) x$fit)
#> $`Regime 1`
#> $`Regime 1`$par
#>           c       gamma 
#> 0.001496696 0.500000000 
#> 
#> $`Regime 1`$objective
#> [1] 3.325991
#> 
#> $`Regime 1`$convergence
#> [1] 0
#> 
#> $`Regime 1`$iterations
#> [1] 60
#> 
#> $`Regime 1`$evaluations
#> function gradient 
#>       81      144 
#> 
#> $`Regime 1`$message
#> [1] "both X-convergence and relative convergence (5)"
#> 
#> 
#> $`Regime 2`
#> $`Regime 2`$par
#>           c       gamma 
#> 0.004074408 0.500000000 
#> 
#> $`Regime 2`$objective
#> [1] 1.106692
#> 
#> $`Regime 2`$convergence
#> [1] 0
#> 
#> $`Regime 2`$iterations
#> [1] 49
#> 
#> $`Regime 2`$evaluations
#> function gradient 
#>       64      121 
#> 
#> $`Regime 2`$message
#> [1] "both X-convergence and relative convergence (5)"
#> 
#> 

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
lapply(fit_temporal[1:2], function(x) x$fit)
#> $`Regime 1`
#> $`Regime 1`$par
#>         a     alpha 
#> 0.4375028 0.3020300 
#> 
#> $`Regime 1`$objective
#> [1] 0.007620233
#> 
#> $`Regime 1`$convergence
#> [1] 0
#> 
#> $`Regime 1`$iterations
#> [1] 14
#> 
#> $`Regime 1`$evaluations
#> function gradient 
#>       20       35 
#> 
#> $`Regime 1`$message
#> [1] "relative convergence (4)"
#> 
#> 
#> $`Regime 2`
#> $`Regime 2`$par
#>         a     alpha 
#> 0.3309062 0.8745089 
#> 
#> $`Regime 2`$objective
#> [1] 0.005785562
#> 
#> $`Regime 2`$convergence
#> [1] 0
#> 
#> $`Regime 2`$iterations
#> [1] 22
#> 
#> $`Regime 2`$evaluations
#> function gradient 
#>       27       52 
#> 
#> $`Regime 2`$message
#> [1] "both X-convergence and relative convergence (5)"
#> 
#> 

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
lapply(fit_sep[1:2], function(x) x$fit)
#> $`Regime 1`
#> $`Regime 1`$par
#>           c       gamma           a       alpha 
#> 0.001539941 0.500000000 0.452841528 0.345950120 
#> 
#> $`Regime 1`$objective
#> [1] 7.62732
#> 
#> $`Regime 1`$convergence
#> [1] 0
#> 
#> $`Regime 1`$iterations
#> [1] 62
#> 
#> $`Regime 1`$evaluations
#> function gradient 
#>       90      290 
#> 
#> $`Regime 1`$message
#> [1] "relative convergence (4)"
#> 
#> 
#> $`Regime 2`
#> $`Regime 2`$par
#>           c       gamma           a       alpha 
#> 0.004107168 0.500000000 0.328388489 0.853896077 
#> 
#> $`Regime 2`$objective
#> [1] 3.302677
#> 
#> $`Regime 2`$convergence
#> [1] 0
#> 
#> $`Regime 2`$iterations
#> [1] 60
#> 
#> $`Regime 2`$evaluations
#> function gradient 
#>       81      281 
#> 
#> $`Regime 2`$message
#> [1] "relative convergence (4)"
#> 
#> 
```
