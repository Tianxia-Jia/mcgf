# Calculate correlation for the general stationary model

Calculate correlation for the general stationary model

## Usage

``` r
..cor_stat(cor_base, lagrangian, lambda, v1, v2, k = 2, h1, h2, u)
```

## Arguments

- cor_base:

  An array of base cross-correlation matrices.

- lagrangian:

  Lagrangian model, `none`, `lagr_tri`, `lagr_askey`, or `lagr_exp`.

- lambda:

  Weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- v1:

  Prevailing wind, u-component.

- v2:

  Prevailing wind, v-component.

- k:

  Scale parameter of \\\\\boldsymbol v\\\\, \\k\>0\\. Default is 2.

- h1:

  Horizontal distance matrix or array.

- h2:

  Vertical distance matrix or array, same dimension as `h1`.

- u:

  Time lag, same dimension as `h1`.

## Value

Correlations for the general stationary model.

## Details

This function is a special case of
[`.cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/dot-cor_stat.md).
It is used inside
[`fit_lagr()`](https://tianxia-jia.github.io/mcgf/reference/fit_lagr.md).
