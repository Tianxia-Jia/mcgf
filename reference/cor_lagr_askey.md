# Calculate Lagrangian correlation of the Askey form

Calculate Lagrangian correlation of the Askey form

## Usage

``` r
cor_lagr_askey(v1, v2, k = 2, h1, h2, u)
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

The Lagrangian correlation function of the Askey form with parameters
\\\boldsymbol v = (v_1, v_2)^\top\in\mathbb{R}^2\\ has the form
\$\$C(\mathbf{h}, u)=\left(1-\dfrac{1}{k\\\boldsymbol v\\}
\left\\\mathbf{h}-u\boldsymbol v\right\\\right)^{3/2}\_+,\$\$ where
\\\\\cdot\\\\ is the Euclidean distance, \\x\_+=\max(x,0), \mathbf{h} =
(\mathrm{h}\_1, \mathrm{h}\_2)^\top\in\mathbb{R}^2\\, and \\k \> 0\\ is
the scale parameter controlling the magnitude of asymmetry in
correlation.

## References

Askey, R. (1973). Radial characteristic functions, Tech. Report No.
1262, Math. Research Center, University of Wisconsin-Madison.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
h1 <- matrix(c(0, -5, 5, 0), nrow = 2)
h2 <- matrix(c(0, -8, 8, 0), nrow = 2)
u <- matrix(0.1, nrow = 2, ncol = 2)
cor_lagr_askey(v1 = 5, v2 = 10, h1 = h1, h2 = h2, u = u)
#>           [,1]      [,2]
#> [1,] 0.9259455 0.4974824
#> [2,] 0.3839919 0.9259455

h1 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
h2 <- array(c(0, -10, 10, 0), dim = c(2, 2, 3))
u <- array(rep(-c(1, 2, 3), each = 4), dim = c(2, 2, 3))
cor_lagr_askey(v1 = 10, v2 = 10, h1 = h1, h2 = h2, u = u)
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.3535534 0.0000000
#> [2,] 1.0000000 0.3535534
#> 
#> , , 2
#> 
#>           [,1] [,2]
#> [1,] 0.0000000    0
#> [2,] 0.3535534    0
#> 
#> , , 3
#> 
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    0    0
#> 
```
