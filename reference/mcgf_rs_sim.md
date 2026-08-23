# Simulate regime-switching Markov chain Gaussian field

Simulate regime-switching Markov chain Gaussian field

## Usage

``` r
mcgf_rs_sim(
  N,
  label,
  base_ls,
  lagrangian_ls,
  par_base_ls,
  par_lagr_ls,
  lambda_ls,
  dists_ls,
  sd_ls,
  lag_ls,
  scale_time = 1,
  init = 0,
  mu_c_ls = list(0),
  mu_p_ls = list(0),
  return_all = FALSE
)
```

## Arguments

- N:

  Sample size.

- label:

  Vector of regime labels of the same length as `N`.

- base_ls:

  List of base model, `sep` or `fs` for now.

- lagrangian_ls:

  List of Lagrangian models. Each element can be `none`, `lagr_tri`,
  `lagr_askey`, or `lagr_exp`.

- par_base_ls:

  List of parameters for the base model.

- par_lagr_ls:

  List of parameters for the Lagrangian model.

- lambda_ls:

  List of weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- dists_ls:

  List of distance matrices or arrays.

- sd_ls:

  List of standard deviation for each location.

- lag_ls:

  List of time lags.

- scale_time:

  Scale of time unit, default is 1. Elements in `lag_ls` are divided by
  `scale_time`.

- init:

  Initial samples, default is 0.

- mu_c_ls, mu_p_ls:

  List of means of current and past.

- return_all:

  Logical; if TRUE the joint covariance matrix, arrays of distances and
  time lag are returned.

## Value

Simulated regime-switching Markov chain Gaussian field with
user-specified covariance structures. The simulation is done by kriging.
The output data is in space-wide format. Each element in `dists_ls` must
contain `h` for symmetric models, and `h1` and `h2` for general
stationary models. `init` can be a scalar or a vector of appropriate
size. List elements in `sd_ls`, `mu_c_ls`, and `mu_p_ls` must be vectors
of appropriate sizes.

## See also

Other simulations of Markov chain Gaussian fields:
[`mcgf_sim()`](https://tianxia-jia.github.io/mcgf/reference/mcgf_sim.md)

## Examples

``` r
par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
par_t <- list(a = 1, alpha = 0.5)
par_base <- list(par_s = par_s, par_t = par_t)
par_lagr <- list(v1 = 5, v2 = 10)
h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
h <- sqrt(h1^2 + h2^2)
dists <- list(h = h, h1 = h1, h2 = h2)

set.seed(123)
label <- sample(1:2, 1000, replace = TRUE)
X <- mcgf_rs_sim(
    N = 1000,
    label = label,
    base_ls = list("sep"),
    lagrangian_ls = list("none", "lagr_tri"),
    lambda_ls = list(0, 0.5),
    par_base_ls = list(par_base),
    par_lagr_ls = list(NULL, par_lagr),
    dists_ls = list(dists, dists)
)
# plot.ts(X[, -1])
```
