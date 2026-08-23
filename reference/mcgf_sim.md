# Simulate Markov chain Gaussian field

Simulate Markov chain Gaussian field

## Usage

``` r
mcgf_sim(
  N,
  base = c("sep", "fs"),
  lagrangian = c("none", "lagr_tri", "lagr_askey", "lagr_exp"),
  par_base,
  par_lagr,
  lambda,
  dists,
  sd = 1,
  lag = 1,
  scale_time = 1,
  horizon = 1,
  init = 0,
  mu_c = 0,
  mu_p = 0,
  return_all = FALSE
)
```

## Arguments

- N:

  Sample size.

- base:

  Base model, `sep` or `fs` for now.

- lagrangian:

  Lagrangian model, `none`, `lagr_tri`, `lagr_askey`, or `lagr_exp`.

- par_base:

  Parameters for the base model (symmetric).

- par_lagr:

  Parameters for the Lagrangian model.

- lambda:

  Weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- dists:

  Distance matrices or arrays.

- sd:

  Standard deviation for each location.

- lag:

  Time lag.

- scale_time:

  Scale of time unit, default is 1. `lag` is divided by `scale_time`.

- horizon:

  Forecast horizon, default is 1.

- init:

  Initial samples, default is 0.

- mu_c, mu_p:

  Means of current and past.

- return_all:

  Logical; if TRUE the joint covariance matrix, arrays of distances and
  time lag are returned.

## Value

Simulated Markov chain Gaussian field with user-specified covariance
structure. The simulation is done by kriging. The output data is in
space-wide format. `dists` must contain `h` for symmetric models, and
`h1` and `h2` for general stationary models. `horizon` controls
forecasting horizon. `sd`, `mu_c`, `mu_p`, and `init` must be vectors of
appropriate sizes.

## See also

Other simulations of Markov chain Gaussian fields:
[`mcgf_rs_sim()`](https://tianxia-jia.github.io/mcgf/reference/mcgf_rs_sim.md)

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
X <- mcgf_sim(
    N = 1000, base = "sep", lagrangian = "lagr_tri", lambda = 0.5,
    par_base = par_base, par_lagr = par_lagr, dists = dists
)
plot.ts(X)

```
