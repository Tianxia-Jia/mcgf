# Calculate exponential correlation

Calculate exponential correlation

## Usage

``` r
cor_exp(x, c, gamma = 1/2, nugget = 0, is.dist = FALSE)
```

## Arguments

- x:

  A numeric vector, matrix, or array.

- c:

  Scale parameter, \\c\>0\\.

- gamma:

  Smooth parameter, \\\gamma\in(0, 1/2\]\\. Default is 1/2.

- nugget:

  The nugget effect \\\in\[0, 1\]\\.

- is.dist:

  Logical; if TRUE, `x` is a distance matrix or an array of distance
  matrices.

## Value

Correlations of the same dimension as `x`.

## Details

The exponential correlation function with scale parameter \\c\\ and
smooth parameter \\\gamma\\ has the form
\$\$C(x)=(1-\text{nugget})\exp(-c\|x\|^{2\gamma})+\text{nugget}\cdot
\delta\_{x=0},\$\$ where \\\delta\_{x=0}\\ is 1 when \\x=0\\ and 0
otherwise.

## References

Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
Geostatistics. Journal of the Royal Statistical Society. Series C
(Applied Statistics), 47(3), 299–350.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
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
cor_exp(x, c = 0.01, gamma = 0.5)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.9512294
#> [2,] 0.9512294 1.0000000

x <- matrix(c(0, 5, 5, 0), nrow = 2)
cor_exp(x, c = 0.01, gamma = 0.5, nugget = 0.3, is.dist = TRUE)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.6658606
#> [2,] 0.6658606 1.0000000
```
