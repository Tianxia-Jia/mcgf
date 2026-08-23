# Calculate correlation for fully symmetric model

Calculate correlation for fully symmetric model

## Usage

``` r
cor_fs(nugget = 0, c, gamma = 1/2, a, alpha, beta = 0, h, u)
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

- beta:

  Interaction parameter, \\\beta\in\[0, 1\]\\.

- h:

  Euclidean distance matrix or array.

- u:

  Time lag, same dimension as `h`.

## Value

Correlations of the same dimension as `h` and `u`.

## Details

The fully symmetric correlation function with interaction parameter
\\\beta\\ has the form \$\$C(\mathbf{h},
u)=\dfrac{1}{(a\|u\|^{2\alpha} + 1)}
\left((1-\text{nugget})\exp\left(\dfrac{-c\\\mathbf{h}\\^{2\gamma}}
{(a\|u\|^{2\alpha}+1)^{\beta\gamma}}\right)+ \text{nugget}\cdot
\delta\_{\mathbf{h}=\boldsymbol 0}\right),\$\$ where \\\\\cdot\\\\ is
the Euclidean distance, and \\\delta\_{x=0}\\ is 1 when \\x=0\\ and 0
otherwise. Here \\\mathbf{h}\in\mathbb{R}^2\\ and \\u\in\mathbb{R}\\. By
default `beta = 0` and it reduces to the separable model.

## References

Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
Space–Time Data, Journal of the American Statistical Association,
97:458, 590-600.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
h <- matrix(c(0, 5, 5, 0), nrow = 2)
u <- matrix(0, nrow = 2, ncol = 2)
cor_fs(c = 0.01, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.9512294
#> [2,] 0.9512294 1.0000000

h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
cor_fs(c = 0.01, gamma = 0.5, a = 1, alpha = 0.5, beta = 0.5, h = h, u = u)
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 1.0000000 0.9512294
#> [2,] 0.9512294 1.0000000
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.5000000 0.4794134
#> [2,] 0.4794134 0.5000000
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.3333333 0.3209070
#> [2,] 0.3209070 0.3333333
#> 
```
