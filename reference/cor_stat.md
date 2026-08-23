# Calculate general stationary correlation.

Calculate general stationary correlation.

## Usage

``` r
cor_stat(
  base = c("sep", "fs"),
  lagrangian = c("none", "lagr_tri", "lagr_askey", "lagr_exp"),
  par_base,
  par_lagr,
  lambda,
  h,
  h1,
  h2,
  u,
  base_fixed = FALSE
)
```

## Arguments

- base:

  Base model, `sep` or `fs` for now. Or correlation matrix/array.

- lagrangian:

  Lagrangian model, `none`, `lagr_tri`, `lagr_askey`, or `lagr_exp`.

- par_base:

  Parameters for the base model (symmetric), used only when
  `base_fixed = FALSE`.

- par_lagr:

  Parameters for the Lagrangian model. Used only when `lagrangian` is
  not `none`.

- lambda:

  Weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- h:

  Euclidean distance matrix or array, used only when
  `base_fixed = FALSE`.

- h1:

  Horizontal distance matrix or array, same dimension as `h`. Used only
  when `lagrangian` is not `none`.

- h2:

  Vertical distance matrix or array, same dimension as `h`. Used only
  when `lagrangian` is not `none`.

- u:

  Time lag, same dimension as `h`.

- base_fixed:

  Logical; if TRUE, `base` is the correlation.

## Value

Correlations for the general stationary model. Same dimension of `base`
if `base_fixed = FALSE`.

## Details

The general station model, a convex combination of a base model and a
Lagrangian model, has the form \$\$C(\mathbf{h},
u)=(1-\lambda)C\_{\text{Base}}(\mathbf{h}, u)+ \lambda
C\_{\text{Lagr}}(\mathbf{h}, u),\$\$ where \\\lambda\\ is the weight of
the Lagrangian term.

If `base_fixed = TRUE`, the correlation is of the form \$\$C(\mathbf{h},
u)=(1-\lambda)C\_{\text{Base}}+ \lambda C\_{\text{Lagr}}(\mathbf{h},
u),\$\$ where `base` is a correlation matrix/array and `par_base` and
`h` are not used.

When `lagrangian = "none"`, `lambda` must be 0.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat_rs()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat_rs.md)

## Examples

``` r
par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
par_t <- list(a = 1, alpha = 0.5)
par_base <- list(par_s = par_s, par_t = par_t)
par_lagr <- list(v1 = 5, v2 = 10)
h1 <- matrix(c(0, 5, -5, 0), nrow = 2)
h2 <- matrix(c(0, 8, -8, 0), nrow = 2)
h <- sqrt(h1^2 + h2^2)
u <- matrix(0.1, nrow = 2, ncol = 2)
cor_stat(
    base = "sep", lagrangian = "lagr_tri", par_base = par_base,
    par_lagr = par_lagr, lambda = 0.8, h = h, h1 = h1, h2 = h2, u = u
)
#>           [,1]      [,2]
#> [1,] 0.9418182 0.5067249
#> [2,] 0.5867249 0.9418182

h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
h <- sqrt(h1^2 + h2^2)
u <- array(rep(c(0.1, 0.2, 0.3), each = 4), dim = c(2, 2, 3))
fit_base <- cor_fs(
    nugget = 0.5, c = 0.01, gamma = 0.5, a = 1, alpha = 0.5,
    beta = 0.0, h = h, u = u
)
par_lagr <- list(v1 = 5, v2 = 10)
cor_stat(
    base = fit_base, lagrangian = "lagr_askey", par_lagr = par_lagr,
    h1 = h1, h2 = h2, u = u, lambda = 0.8, base_fixed = TRUE
)
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.9225746 0.3899184
#> [2,] 0.4807108 0.9225746
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.8497186 0.3405980
#> [2,] 0.5219630 0.8497186
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.7807752 0.2944717
#> [2,] 0.5659495 0.7807752
#> 
```
