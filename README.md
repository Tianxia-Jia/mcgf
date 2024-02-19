
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mcgf

<!-- badges: start -->

[![R-CMD-check](https://github.com/tianxia-jia/mcgf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tianxia-jia/mcgf/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/mcgf)](https://CRAN.R-project.org/package=mcgf)
<!-- badges: end -->

The goal of `mcgf` is to provide easy-to-use functions for simulating
and fitting covariance models. It provides functions for simulating
(regime-switching) Markov chain Gaussian fields with covariance
functions of the Gneiting class by simple kriging. Parameter estimation
methods such as weighted least squares and maximum likelihood estimation
are available. Below is an example of simulating and estimation
parameters for an MCGF.

## Installation

You can install the development version of mcgf from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tianxia-jia/mcgf")
```

## Data Simulation

To simulate an MCGF with fully symmetric covariance structure, we begin
with simulating 10 locations randomly.

``` r
library(mcgf)
set.seed(123)
h <- rdists(10)
```

Next, we simulate an MCGF with the general stationary covariance
structure. In this example the covariance structure is a convex
combination of a base separable model and a Lagrangian model account for
asymmetry.

``` r
N <- 1000
lag <- 5

par_base <- list(
    par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.8)
)
par_lagr <- list(v1 = 200, v2 = 200, k = 2)

sim1 <- mcgf_sim(
    N = N,
    base = "sep",
    lagrangian = "lagr_tri",
    par_base = par_base,
    par_lagr = par_lagr,
    lambda = 0.2,
    dists = h,
    lag = lag
)
sim1 <- sim1[-c(1:(lag + 1)), ]
rownames(sim1) <- 1:nrow(sim1)

sim1 <- list(data = sim1, dists = h)
```

## Parameter Estimation

### Create an `mcgf` object

To estimate parameters, we need to calculate auto-correlations and
cross-correlations. Let’s first create an `mcgf` object. The `mcgf`
class extends the `data.frame` with more attributes.

``` r
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
```

Then the acfs and ccfs can be added to this object as follows.

``` r
sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = lag)
sim1_mcgf <- add_ccfs(sim1_mcgf, lag_max = lag)
```

### Estimate base model

To perform parameter estimation, we can start with estimating the
parameters for spatial and temporal models.

``` r
fit_spatial <- fit_base(
    sim1_mcgf,
    model = "spatial",
    lag = lag,
    par_init = c(c = 0.001, gamma = 0.5),
    par_fixed = c(nugget = 0)
)
fit_spatial$fit
#> $par
#>           c       gamma 
#> 0.001160802 0.500000000 
#> 
#> $objective
#> [1] 1.640593
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 8
#> 
#> $evaluations
#> function gradient 
#>       21       20 
#> 
#> $message
#> [1] "both X-convergence and relative convergence (5)"
```

``` r
fit_temporal <- fit_base(
    sim1_mcgf,
    model = "temporal",
    lag = lag,
    par_init = c(a = 0.3, alpha = 0.5)
)
fit_temporal$fit
#> $par
#>         a     alpha 
#> 0.6528906 0.7560970 
#> 
#> $objective
#> [1] 0.004306706
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 18
#> 
#> $evaluations
#> function gradient 
#>       23       43 
#> 
#> $message
#> [1] "relative convergence (4)"
```

Alternatively, we can fit the separable model all at once:

``` r
fit_sep <- fit_base(
    sim1_mcgf,
    model = "sep",
    lag = lag,
    par_init = c(
        c = 0.001,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.5
    ),
    par_fixed = c(nugget = 0)
)
fit_sep$fit
#> $par
#>           c       gamma           a       alpha 
#> 0.001154864 0.500000000 0.624551338 0.735490605 
#> 
#> $objective
#> [1] 3.488305
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 18
#> 
#> $evaluations
#> function gradient 
#>       49       88 
#> 
#> $message
#> [1] "relative convergence (4)"
```

we can also estimate the parameters using MLE:

``` r
fit_sep2 <- fit_base(
    sim1_mcgf,
    model = "sep",
    lag = lag,
    par_init = c(
        c = 0.001,
        gamma = 0.5,
        a = 0.5,
        alpha = 0.5
    ),
    par_fixed = c(nugget = 0),
    method = "mle",
)
fit_sep2$fit
#> $par
#>           c       gamma           a       alpha 
#> 0.001197799 0.500000000 0.804621207 1.000000000 
#> 
#> $objective
#> [1] -11520.04
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 17
#> 
#> $evaluations
#> function gradient 
#>       55       78 
#> 
#> $message
#> [1] "relative convergence (4)"
```

Now we will add the base model to `x_mcgf`:

``` r
sim1_mcgf <- add_base(sim1_mcgf, fit_base = fit_sep)
```

To print the current model, we do

``` r
model(sim1_mcgf)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma           a       alpha      nugget 
#> 0.001154864 0.500000000 0.624551338 0.735490605 0.000000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> ----------------------------------------
#>               Lagrangian
#> ----------------------------------------
#> - Lagrangian model: 
#> - Parameters:
#> NULL
#> 
#> - Fixed parameters:
#> NULL
#> 
#> - Parameter estimation method: 
#> 
#> - Optimization function:
```

### Estimate the Lagrangian model

Similarly, we can estimate the parameters for the Lagrangian component
by

``` r
fit_lagr <- fit_lagr(
    sim1_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 300, v2 = 300, lambda = 0.15),
    par_fixed = c(k = 2)
)
fit_lagr$fit
#> $par
#>      lambda          v1          v2 
#>   0.1757035 232.0852117 203.8869305 
#> 
#> $objective
#> [1] 1.627017
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 32
#> 
#> $evaluations
#> function gradient 
#>       35      126 
#> 
#> $message
#> [1] "relative convergence (4)"
```

We can add the Lagrangian model by

``` r
sim1_mcgf <- add_lagr(sim1_mcgf, fit_lagr = fit_lagr)
```

Finally we may print the final model:

``` r
model(sim1_mcgf)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Base model: sep 
#> - Parameters:
#>           c       gamma           a       alpha      nugget 
#> 0.001154864 0.500000000 0.624551338 0.735490605 0.000000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> ----------------------------------------
#>               Lagrangian
#> ----------------------------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>      lambda          v1          v2           k 
#>   0.1757035 232.0852117 203.8869305   2.0000000 
#> 
#> - Fixed parameters:
#> k 
#> 2 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb
```

### Kriging forecast

This package provides kriging forecasts (and intervals) for empirical,
base, and general stationary models.

``` r
# Empirical model
fit_emp <-
    krige(sim1_mcgf,
        model = "empirical",
        interval = TRUE
    )
rmse_emp <- sqrt(mean(colMeans((sim1_mcgf - fit_emp$fit)^2, na.rm = T)))

# Base separable model
fit_base <-
    krige(sim1_mcgf,
        model = "base",
        interval = TRUE
    )
rmse_base <-
    sqrt(mean(colMeans((sim1_mcgf - fit_base$fit)^2, na.rm = T)))

# Stationary model
fit_stat <-
    krige(sim1_mcgf,
        model = "all",
        interval = TRUE
    )
rmse_stat <-
    sqrt(mean(colMeans((sim1_mcgf - fit_stat$fit)^2, na.rm = T)))

rmse <- c(rmse_emp, rmse_base, rmse_stat)
names(rmse) <- c("Empirical", "Separable", "Stationary")
rmse
#>  Empirical  Separable Stationary 
#>  0.7212175  0.7685016  0.7355458
```
