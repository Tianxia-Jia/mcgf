# Parameter estimation for Lagrangian correlation functions for an `mcgf_rs` object.

Parameter estimation for Lagrangian correlation functions for an
`mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf_rs'
fit_lagr(
  x,
  model_ls,
  method_ls = "wls",
  optim_fn_ls = "nlminb",
  par_fixed_ls = list(NULL),
  par_init_ls,
  lower_ls = list(NULL),
  upper_ls = list(NULL),
  other_optim_fn_ls = list(NULL),
  dists_base_ls,
  dists_lagr_ls = list(NULL),
  rs = TRUE,
  ...
)
```

## Arguments

- x:

  An `mcgf_rs` object containing attributes `dists`, `acfs`, `ccfs`, and
  `sds`. `x` must have been passed to
  [`add_base()`](https://tianxia-jia.github.io/mcgf/reference/add_base.md)
  or `base<-`

- model_ls:

  List of base models, each element must be one of `lagr_tri`,
  `lagr_askey`, `lagr_exp`, or `none`. If `none`, NULLs are returned.

- method_ls:

  List of parameter estimation methods, weighted least square (`wls`) or
  maximum likelihood estimation (`mle`).

- optim_fn_ls:

  List of optimization functions, each element must be one of `nlminb`,
  `optim`, `other`. When use `other`, supply `other_optim_fn_ls`

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

- dists_lagr_ls:

  List of distance matrices/arrays. Used when `dists_base` is FALSE. If
  NULL, `dists(x)` is used.

- rs:

  Logical; if TRUE `x` is treated as a regime-switching, FALSE if the
  parameters need to be estimated in a non-regime-switching setting.

- ...:

  Additional arguments passed to `optim_fn`.

## Value

A list containing outputs from optimization functions of `optim_fn`.

## Details

This functions is the regime-switching variant of
[`fit_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf.md).
Arguments are in lists. The length of arguments that end in `_ls` must
be 1 or the same as the number of regimes in `x`. If the length of an
argument is 1, then it is set the same for all regimes. Refer to
[`fit_lagr.mcgf()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.mcgf.md)
for more details of the arguments.

Note that both `wls` and `mle` are heuristic approaches when `x`
contains observations from a subset of the discrete spatial domain,
though estimation results are close to that using the full spatial
domain for large sample sizes.

Since parameters for the base model and the Lagrangian model are
estimated sequentially, more accurate estimation may be obtained if the
full model is fitted all at once.

## See also

Other functions on fitting an mcgf_rs:
[`add_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_base.mcgf_rs.md),
[`add_lagr.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/add_lagr.mcgf_rs.md),
[`fit_base.mcgf_rs()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.mcgf_rs.md),
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
lapply(fit_lagr_rs[1:2], function(x) x$fit)
#> $`Regime 1`
#> $`Regime 1`$par
#>        v1        v2 
#> -106.5181  119.4434 
#> 
#> $`Regime 1`$objective
#> [1] 5.137749
#> 
#> $`Regime 1`$convergence
#> [1] 0
#> 
#> $`Regime 1`$iterations
#> [1] 26
#> 
#> $`Regime 1`$evaluations
#> function gradient 
#>       29       61 
#> 
#> $`Regime 1`$message
#> [1] "relative convergence (4)"
#> 
#> 
#> $`Regime 2`
#> $`Regime 2`$par
#>       v1       v2 
#> 210.1156 242.3082 
#> 
#> $`Regime 2`$objective
#> [1] 5.902847
#> 
#> $`Regime 2`$convergence
#> [1] 0
#> 
#> $`Regime 2`$iterations
#> [1] 29
#> 
#> $`Regime 2`$evaluations
#> function gradient 
#>       31       77 
#> 
#> $`Regime 2`$message
#> [1] "relative convergence (4)"
#> 
#> 
```
