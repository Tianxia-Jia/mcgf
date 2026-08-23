# Compute the objective for wls method

Compute the objective for wls method

## Usage

``` r
obj_wls(par, cor_fn, cor_emp, par_fixed)
```

## Arguments

- par:

  Parameters of `cor_fn`.

- cor_fn:

  Correlation function.

- cor_emp:

  Empirical correlations.

- par_fixed:

  Fixed parameters of `cor_fn`.

## Value

The objective of weighted least squares.
