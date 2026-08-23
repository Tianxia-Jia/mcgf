# Compute the objective for mle method

Compute the objective for mle method

## Usage

``` r
obj_mle(par, cor_fn, x, lag, par_fixed)
```

## Arguments

- par:

  Parameters of `cor_fn`.

- cor_fn:

  Correlation function

- x:

  An `mcgf` or `mcgf_rs` object

- lag:

  Time lag.

- par_fixed:

  Fixed parameters of `cor_fn`.

## Value

A numeric value proportional to the negative approximate conditional
log-likelihood, to be minimized during parameter estimation.
