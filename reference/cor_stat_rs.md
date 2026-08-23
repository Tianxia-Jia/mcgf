# Calculate general stationary correlation.

Calculate general stationary correlation.

## Usage

``` r
cor_stat_rs(
  n_regime,
  base_ls,
  lagrangian_ls,
  par_base_ls,
  par_lagr_ls,
  lambda_ls,
  h_ls,
  h1_ls,
  h2_ls,
  u_ls,
  base_fixed = FALSE
)
```

## Arguments

- n_regime:

  Integer, number of regimes.

- base_ls:

  List of base model, `sep` or `fs` for now. Or list of correlation
  matrices/arrays.

- lagrangian_ls:

  List of Lagrangian models. Each element can be `none`, `lagr_tri`,
  `lagr_askey`, or `lagr_exp`.

- par_base_ls:

  List of parameters for the base model, used only when
  `base_fixed = FALSE`.

- par_lagr_ls:

  List of parameters for the Lagrangian model. Used only when
  `lagrangian_ls` is not `none`.

- lambda_ls:

  List of weight of the Lagrangian term, \\\lambda\in\[0, 1\]\\.

- h_ls:

  List of Euclidean distance matrix or array, used only when
  `base_fixed = FALSE`.

- h1_ls:

  List of horizontal distance matrix or array, same dimension as `h_ls`.
  Used only when `lagrangian_ls` is not `none`.

- h2_ls:

  List of vertical distance matrix or array, same dimension as `h_ls`.
  Used only when `lagrangian_ls` is not `none`.

- u_ls:

  List of time lag, same dimension as `h_ls`.

- base_fixed:

  Logical; if TRUE, `base_ls` is the list of correlation.

## Value

Correlations for the general stationary model. Same dimension of
`base_ls` if `base_fixed = TRUE`.

## Details

It gives a list of general stationary correlation for `n_regime`
regimes. See
[cor_stat](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md) for
the model details. Model parameters are lists of length 1 or `n_regime`.
When length is 1, same values are used for all regimes. If
`base_fixed = TRUE`, the base is a list of correlation and `par_base_ls`
and `h_ls` are not used.

## See also

Other correlation functions:
[`cor_cauchy()`](https://tianxia-jia.github.io/mcgf/reference/cor_cauchy.md),
[`cor_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_exp.md),
[`cor_fs()`](https://tianxia-jia.github.io/mcgf/reference/cor_fs.md),
[`cor_lagr_askey()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_askey.md),
[`cor_lagr_exp()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_exp.md),
[`cor_lagr_tri()`](https://tianxia-jia.github.io/mcgf/reference/cor_lagr_tri.md),
[`cor_sep()`](https://tianxia-jia.github.io/mcgf/reference/cor_sep.md),
[`cor_stat()`](https://tianxia-jia.github.io/mcgf/reference/cor_stat.md)

## Examples

``` r
# Fit general stationary model with sep base.
par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5)
par_t <- list(a = 1, alpha = 0.5)
par_base <- list(par_s = par_s, par_t = par_t)
h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
h <- sqrt(h1^2 + h2^2)
u <- array(rep(c(1, 2, 3), each = 4), dim = c(2, 2, 3))
cor_stat_rs(
    n_regime = 2,
    base_ls = list("sep"),
    lagrangian_ls = list("none", "lagr_tri"),
    par_base_ls = list(par_base),
    par_lagr_ls = list(NULL, list(v1 = 10, v2 = 20)),
    lambda_ls = list(0, 0.2),
    h_ls = list(h),
    h1_ls = list(NULL, h1),
    h2_ls = list(NULL, h2),
    u_ls = list(u, u + 1)
)
#> [[1]]
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.5000000 0.2274934
#> [2,] 0.2274934 0.5000000
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.3333333 0.1516622
#> [2,] 0.1516622 0.3333333
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.2500000 0.1137467
#> [2,] 0.1137467 0.2500000
#> 
#> 
#> [[2]]
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.2666667 0.1213298
#> [2,] 0.1633298 0.2666667
#> 
#> , , 2
#> 
#>            [,1]       [,2]
#> [1,] 0.20000000 0.09099735
#> [2,] 0.09099735 0.20000000
#> 
#> , , 3
#> 
#>            [,1]       [,2]
#> [1,] 0.16000000 0.07279788
#> [2,] 0.07279788 0.16000000
#> 
#> 

# Fit general stationary model given fs as the base model.
h1 <- array(c(0, 5, -5, 0), dim = c(2, 2, 3))
h2 <- array(c(0, 8, -8, 0), dim = c(2, 2, 3))
h <- sqrt(h1^2 + h2^2)
u <- array(rep(c(0.1, 0.2, 0.3), each = 4), dim = c(2, 2, 3))
fit_base <- cor_fs(
    nugget = 0.5, c = 0.01, gamma = 0.5, a = 1, alpha = 0.5,
    beta = 0.0, h = h, u = u
)
par_lagr <- list(v1 = 5, v2 = 10)
cor_stat_rs(
    n_regime = 2,
    par_lagr_ls = list(par_lagr),
    h1_ls = list(h1),
    h2_ls = list(h2),
    u_ls = list(u, u + 1),
    lambda_ls = list(0, 0.8),
    base_ls = list(fit_base),
    lagrangian_ls = list("lagr_tri", "lagr_askey"),
    base_fixed = TRUE
)
#> [[1]]
#> , , 1
#> 
#>           [,1]      [,2]
#> [1,] 0.9090909 0.4136243
#> [2,] 0.4136243 0.9090909
#> 
#> , , 2
#> 
#>           [,1]      [,2]
#> [1,] 0.8333333 0.3791556
#> [2,] 0.3791556 0.8333333
#> 
#> , , 3
#> 
#>           [,1]      [,2]
#> [1,] 0.7692308 0.3499898
#> [2,] 0.3499898 0.7692308
#> 
#> 
#> [[2]]
#> , , 1
#> 
#>           [,1]       [,2]
#> [1,] 0.4233135 0.08671162
#> [2,] 0.7251899 0.42331352
#> 
#> , , 2
#> 
#>           [,1]       [,2]
#> [1,] 0.3690524 0.07583112
#> [2,] 0.6650994 0.36905244
#> 
#> , , 3
#> 
#>           [,1]       [,2]
#> [1,] 0.3194964 0.06999796
#> [2,] 0.6069045 0.31949639
#> 
#> 
```
