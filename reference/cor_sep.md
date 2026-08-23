# Calculate correlation for separable model

Calculate correlation for separable model

## Usage

``` r
cor_sep(
  spatial = c("exp", "cauchy"),
  temporal = c("exp", "cauchy"),
  par_s,
  par_t,
  h,
  u
)
```

## Arguments

- spatial:

  Pure spatial model, `exp` or `cauchy` for now.

- temporal:

  Pure temporal model, `exp` or `cauchy` for now.

- par_s:

  Parameters for the pure spatial model. Nugget effect supported.

- par_t:

  Parameters for the pure temporal model.

- h:

  Euclidean distance matrix or array.

- u:

  Time lag, same dimension as `h`.

## Value

Correlations of the same dimension as `h` and `u`.

## Details

The separable model is the product of a pure temporal model, \\C_T(u)\\,
and a pure spatial model, \\C_S(\mathbf{h})\\. It is of the form
\$\$C(\mathbf{h}, u)=C\_{T}(u)
\left\[(1-\text{nugget})C\_{S}(\mathbf{h})+\text{nugget}
\delta\_{\mathbf{h}=0}\right\],\$\$ where \\\delta\_{x=0}\\ is 1 when
\\x=0\\ and 0 otherwise. Here \\\mathbf{h}\in\mathbb{R}^2\\ and
\\u\in\mathbb{R}\\. Now only exponential and Cauchy correlation models
are available.

## References

Gneiting, T. (2002). Nonseparable, Stationary Covariance Functions for
Space–Time Data, Journal of the American Statistical Association,
97:458, 590-600.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
h <- matrix(c(0, 5, 5, 0), nrow = 2)
par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
u <- matrix(0, nrow = 2, ncol = 2)
par_t <- list(a = 1, alpha = 0.5)
cor_sep(
    spatial = "exp", temporal = "cauchy", par_s = par_s, par_t = par_t,
    h = h, u = u
)
#>           [,1]      [,2]
#> [1,] 1.0000000 0.4756147
#> [2,] 0.4756147 1.0000000

h <- array(c(0, 5, 5, 0), dim = c(2, 2, 3))
par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
u <- array(rep(0:2, each = 4), dim = c(2, 2, 3))
par_t <- list(a = 1, alpha = 0.5)
cor_sep(
    spatial = "exp", temporal = "cauchy", par_s = par_s, par_t = par_t,
    h = h, u = u
)
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 1.0000000 0.4756147
#> [2,] 0.4756147 1.0000000
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.5000000 0.2378074
#> [2,] 0.2378074 0.5000000
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.3333333 0.1585382
#> [2,] 0.1585382 0.3333333
#> 
```
