# Calculate Lagrangian correlation of the exponential form

Calculate Lagrangian correlation of the exponential form

## Usage

``` r
cor_lagr_exp(v1, v2, k = 2, h1, h2, u)
```

## Arguments

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

Correlations of the same dimension as `h1`.

## Details

The Lagrangian correlation function of the exponential form with
parameters \\\boldsymbol v = (v_1, v_2)^\top\in\mathbb{R}^2\\ has the
form \$\$C(\mathbf{h}, u)=\exp\left(-\dfrac{1}{k\\\boldsymbol v\\}
\left\\\mathbf{h}-u\boldsymbol v\right\\\right),\$\$ where \\\\\cdot\\\\
is the Euclidean distance, \\\mathbf{h} = (\mathrm{h}\_1,
\mathrm{h}\_2)^\top\in\mathbb{R}^2\\, and \\k \> 0\\ is the scale
parameter controlling the magnitude of asymmetry in correlation.

## References

Diggle, P. J., Tawn, J. A., & Moyeed, R. A. (1998). Model-Based
Geostatistics. Journal of the Royal Statistical Society. Series C
(Applied Statistics), 47(3), 299–350.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
h1 <- matrix(c(0, -5, 5, 0), nrow = 2)
h2 <- matrix(c(0, -8, 8, 0), nrow = 2)
u <- matrix(0.1, nrow = 2, ncol = 2)
cor_lagr_exp(v1 = 5, v2 = 10, h1 = h1, h2 = h2, u = u)
#>           [,1]      [,2]
#> [1,] 0.9512294 0.6892468
#> [2,] 0.6239413 0.9512294

h1 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
h2 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
u <- array(rep(-c(1, 2, 3), each = 4), dim = c(2, 2, 3))
cor_lagr_exp(v1 = 10, v2 = 10, h1 = h1, h2 = h2, u = u)
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.6065307 0.3678794
#> [2,] 1.0000000 0.6065307
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.3678794 0.2231302
#> [2,] 0.6065307 0.3678794
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.2231302 0.1353353
#> [2,] 0.3678794 0.2231302
#> 
```
