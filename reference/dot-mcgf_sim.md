# Simulate Markov chain Gaussian field

Simulate Markov chain Gaussian field

## Usage

``` r
.mcgf_sim(
  N,
  base,
  lagrangian,
  par_base,
  par_lagr,
  lambda,
  dists,
  sd,
  lag,
  scale_time = 1,
  horizon = 1,
  init = 0,
  mu_c,
  mu_p,
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
