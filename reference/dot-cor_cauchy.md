# Calculate Cauchy correlation

Calculate Cauchy correlation

## Usage

``` r
.cor_cauchy(x, a, alpha, nu = 1, nugget = 0)
```

## Arguments

- x:

  A numeric vector, matrix, or array.

- a:

  Scale parameter, \\a\>0\\.

- alpha:

  Smooth parameter, \\\alpha\in(0, 1\]\\.

- nu:

  Power parameter, \\\nu\>0\\. Default is 1.

- nugget:

  The nugget effect \\\in\[0, 1\]\\.

## Value

Correlations of the same dimension as `x`.

## Details

The Cauchy correlation function with scale parameter \\a\\ and smooth
parameter \\\alpha\\ has the form
\$\$C(x)=(1-\text{nugget})(a\|x\|^{2\alpha} +
1)^{-\nu}+\text{nugget}\cdot \delta\_{x=0},\$\$ where \\\delta\_{x=0}\\
is 1 when \\x=0\\ and 0 otherwise.

## References

Gneiting, T., and Schlather, M. (2004). Stochastic Models That Separate
Fractal Dimension and the Hurst Effect. SIAM Review, 46(2), 269–282.
