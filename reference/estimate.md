# Optimization for wls method

Optimization for wls method

## Usage

``` r
estimate(par_init, method, optim_fn, cor_fn, par_fixed, lower, upper, ...)
```

## Arguments

- par_init:

  Initial values for parameters to be optimized.

- method:

  Parameter estimation method. "wls" or "mle".

- optim_fn:

  Optimization function.

- cor_fn:

  Correlation function.

- par_fixed:

  Fixed parameters of `cor_fn`.

- lower:

  Lower bound.

- upper:

  Upper bound.

- ...:

  Additional arguments passed to `optim_fn`, `obj_wls` or `obj_mle`.

## Value

A list outputted from optimization functions of `optim_fn`.
