
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcgf

<!-- badges: start -->
<!-- badges: end -->

The goal of `mcgf` is to provide easy-to-use functions for simulating
and fitting covariance models. This package has functions to simulate
samples for Markov chain Gaussian fields with covariance functions of
the Gneiting class by simple kriging. Parameter estimation methods such
as weighted least squares and maximum likelihood estimation are
available.

## Installation

You can install the development version of mcgf from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tianxia-jia/mcgf")
```

## Data Simulation

To simulate a MCGF with fully symmetric covariance structure, we first
create distance arrays for 10 locations:

``` r
library(mcgf)
set.seed(123)
dists <- rdists(10)
#> The legacy packages maptools, rgdal, and rgeos, underpinning the sp package,
#> which was just loaded, will retire in October 2023.
#> Please refer to R-spatial evolution reports for details, especially
#> https://r-spatial.org/r/2023/05/15/evolution4.html.
#> It may be desirable to make the sf package available;
#> package maintainers should consider adding sf to Suggests:.
#> The sp package is now running under evolution status 2
#>      (status 2 uses the sf package in place of rgdal)
```

Next, we simulate a MCGF.

``` r
library(mcgf)
par_s <- list(nugget = 0.2, c = 0.001, gamma = 0.3)
par_t <- list(a = 5, alpha = 0.5)
par_base <- list(par_s = par_s, par_t = par_t)

set.seed(123)
X <- mcgf_sim(N = 1000, base = "sep", lagrangian = "none", par_base = par_base,
              dists = dists, lag = 10)
plot.ts(X)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />

## Parameter Estimation

To estimate parameters using weighted least squares method, we need to
calculate autocorrelations and cross-autocorrelations. Let’s first
create a `mcgf` object. The `mcgf` class extends the `data.frame` with
more attributes.

``` r
x_mcgf <- mcgf(X, dists = dists)
#> `time` not provided, assuming rows are equally spaced temporally.
```

Then the acfs and ccfs can be added to this object as follows.

``` r
x_mcgf <- add_acfs(x = x_mcgf, lag_max = 10)
x_mcgf <- add_ccfs(x = x_mcgf, lag_max = 10)
```

To find the weighted least squares, we can run

``` r
fit_base(x = x_mcgf, lag = 10, model = "sep", optim_fn = "nlminb",
         par_init = list(a = 1, alpha = 0.5, nugget = 0.1, c = 0.001, 
                         gamma = 0.5), method = "wls")
#> $par
#>            c        gamma       nugget            a        alpha 
#> 3.226771e-05 5.000000e-01 2.056779e-01 4.068848e+00 5.096764e-01 
#> 
#> $objective
#> [1] 1.000885
#> 
#> $convergence
#> [1] 1
#> 
#> $iterations
#> [1] 150
#> 
#> $evaluations
#> function gradient 
#>      192      821 
#> 
#> $message
#> [1] "iteration limit reached without convergence (10)"
#> 
#> $par_names
#> [1] "c"      "gamma"  "nugget" "a"      "alpha"
```

To find the maximum likelihood estimates, we run

``` r
fit_base(x = x_mcgf, lag = 10, model = "sep", optim_fn = "nlminb",
         par_init = list(a = 1, alpha = 0.5, nugget = 0.1, c = 0.001, 
                         gamma = 0.5), method = "mle")
#> $par
#>           c       gamma      nugget           a       alpha 
#> 0.008520959 0.167558413 0.172344741 4.945943574 0.510839525 
#> 
#> $objective
#> [1] -1927.716
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 76
#> 
#> $evaluations
#> function gradient 
#>      106      527 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
#> $par_names
#> [1] "c"      "gamma"  "nugget" "a"      "alpha"
```

<!-- ## Regime-Switching MCGF -->
<!-- Regime-switching MCGF can be simulated as -->
<!-- ```{r, warning=FALSE,error=FALSE} -->
<!-- par_s <- list(nugget = 0.5, c = 0.01, gamma = 0.5) -->
<!-- par_t <- list(a = 1, alpha = 0.5) -->
<!-- par_base <- list(par_s = par_s, par_t = par_t) -->
<!-- par_lagr <- list(v1 = 5, v2 = 10) -->
<!-- set.seed(123) -->
<!-- labels <- sample(1:2, 1000, replace = TRUE) -->
<!-- X <- mcgf_rs_sim(N = 1000, -->
<!--                  labels = labels, -->
<!--                  base_ls = list("sep"), -->
<!--                  lagrangian_ls = list("none", "none"), -->
<!--                  lambda_ls = list(0, 0.5), -->
<!--                  par_base_ls = list(par_base), -->
<!--                  par_lagr_ls = list(NULL, par_lagr), -->
<!--                  dists_ls = list(dists, dists)) -->
<!-- plot.ts(X[, -1]) -->
<!-- ``` -->
