# Calculate correlation for separable model

Calculate correlation for separable model

## Usage

``` r
..cor_sep(nugget, c, gamma = 1/2, a, alpha, h, u)
```

## Arguments

- nugget:

  The nugget effect \\\in\[0, 1\]\\.

- c:

  Scale parameter of `cor_exp`, \\c\>0\\.

- gamma:

  Smooth parameter of `cor_exp`, \\\gamma\in(0, 1/2\]\\.

- a:

  Scale parameter of `cor_cauchy`, \\a\>0\\.

- alpha:

  Smooth parameter of `cor_cauchy`, \\\alpha\in(0, 1\]\\.

- h:

  Euclidean distance matrix or array.

- u:

  Time lag, same dimension as `h`.

## Value

Correlations for separable model.

## Details

This function is a special case of
[`.cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/dot-cor_fs.md).
It is used inside
[`fit_base()`](https://tianxia-jia.github.io/mcgf/reference/fit_base.md).
