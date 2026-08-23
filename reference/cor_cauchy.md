# Calculate Cauchy correlation

Calculate Cauchy correlation

## Usage

``` r
cor_cauchy(x, a, alpha, nu = 1, nugget = 0, is.dist = FALSE)
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

- is.dist:

  Logical; if TRUE, `x` is a distance matrix or an array of distance
  matrices.

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

## See also

Other correlation functions:
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
x <- matrix(c(0, 5, 5, 0), nrow = 2)
cor_cauchy(x, a = 1, alpha = 0.5)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.1666667
#> [2,] 0.1666667 1.0000000

x <- matrix(c(0, 5, 5, 0), nrow = 2)
cor_cauchy(x, a = 1, alpha = 0.5, nugget = 0.3, is.dist = TRUE)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.1166667
#> [2,] 0.1166667 1.0000000
```
