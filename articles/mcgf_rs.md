# Steps for fitting an \`mcgf_rs\` object

## Introduction

The `mcgf` package contains useful functions to simulate Markov chain
Gaussian fields (MCGF) and regime-switching Markov chain Gaussian fields
(RS-MCGF) with covariance structures of the Gneiting class (Gneiting
2002), including regime-switching spatio-temporal covariance models
developed for short-term wind forecasting (Jia and Sezer 2025). It also
provides useful tools to estimate the parameters by weighted least
squares (WLS) and maximum likelihood estimation (MLE). The `mcgf_rs`
function can be used to fit covariance models and obtain Kriging
forecasts. A typical workflow for fitting a regime-switching MCGF is
given below.

1.  Create an `mcgf_rs` object by providing a dataset and the
    corresponding locations/distances
2.  Calculate auto- and cross-correlations.
3.  Fit the base covariance model which is a fully symmetric model.
4.  Fit the Lagrangian model to account for asymmetry, if necessary.
5.  Obtain Kriging forecasts.
6.  Obtain Kriging forecasts for new locations given their coordinates.

We will demonstrate the use of `mcgf_rs` by an example below.

## Simulated Example

Here we will generate a regime-switching MCGF process of size 5,000 with
the general stationary covariance structure. More specifically, we
suppose there are two regimes each with different prevailing winds but
the same fully symmetric model. We will begin with simulating the regime
process first.

### Regime process

We start with simulating a transition matrix for a two-state Markov
chain with high probabilities of staying in the current regime. The
simulated transition matrix and regime process are given below.

``` r

library(mcgf)
K <- 2
N <- 5000
lag <- 5

set.seed(123)
tran_mat <- matrix(rnorm(K^2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
tran_mat
#>            [,1]       [,2]
#> [1,] 0.92650859 0.07349141
#> [2,] 0.04934827 0.95065173

regime <- rep(NA, N)
regime[1] <- 1

for (n in 2:(N)) {
    regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
}
table(regime)
#> regime
#>    1    2 
#> 1993 3007
```

### Data Simulation

#### Locations

We generate the coordinates of 25 locations below and treat the last 5
locations as unobserved locations.

``` r

set.seed(123)
x <- stats::rnorm(25, -110)
y <- stats::rnorm(25, 50)
locations_all <- cbind(lon = x, lat = y)
locations <- locations_all[1:20, ]
locations_new <- locations_all[-c(1:20), ]
locations_all
#>             lon      lat
#>  [1,] -110.5605 48.31331
#>  [2,] -110.2302 50.83779
#>  [3,] -108.4413 50.15337
#>  [4,] -109.9295 48.86186
#>  [5,] -109.8707 51.25381
#>  [6,] -108.2849 50.42646
#>  [7,] -109.5391 49.70493
#>  [8,] -111.2651 50.89513
#>  [9,] -110.6869 50.87813
#> [10,] -110.4457 50.82158
#> [11,] -108.7759 50.68864
#> [12,] -109.6402 50.55392
#> [13,] -109.5992 49.93809
#> [14,] -109.8893 49.69404
#> [15,] -110.5558 49.61953
#> [16,] -108.2131 49.30529
#> [17,] -109.5021 49.79208
#> [18,] -111.9666 48.73460
#> [19,] -109.2986 52.16896
#> [20,] -110.4728 51.20796
#> [21,] -111.0678 48.87689
#> [22,] -110.2180 49.59712
#> [23,] -111.0260 49.53334
#> [24,] -110.7289 50.77997
#> [25,] -110.6250 49.91663
```

#### Data simulation

To simulate the 5th-order RS-MCGF process, we calculate the distance
matrices for all locations first and then use `mcgf_rs_sim` to simulate
such process.

``` r

# simulate RS MCGF
par_base <- list(
    par_s = list(nugget = 0, c = 0.005, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)
par_lagr1 <- list(v1 = 100, v2 = 100, k = 2)
par_lagr2 <- list(v1 = 50, v2 = 100, k = 2)

h <- find_dists_new(locations, locations_new)

set.seed(123)
data_all <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"),
    lagrangian_ls = list("lagr_tri"),
    par_base_ls = list(par_base),
    par_lagr_ls = list(par_lagr1, par_lagr2),
    lambda_ls = list(0.3, 0.3),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
data_all <- data_all[-c(1:(lag + 1)), ]
rownames(data_all) <- 1:nrow(data_all)
head(data_all)
#>   regime          1           2          3           4           5          6
#> 1      1 -0.4131534 -0.21405685  0.7229975 -0.08322733  0.07023563  1.1604699
#> 2      1 -1.7574429 -0.06411826  0.3964733 -1.75421279  0.80900801  0.9012895
#> 3      1 -0.8264216 -1.02490415 -0.4809606 -0.52653010 -0.49018553  0.4388830
#> 4      1  0.4089212 -0.72910670 -1.6972114 -0.01897436 -1.38506809 -1.1222586
#> 5      1 -0.4281591  0.03573133 -0.4868917 -0.58029926 -0.71498960 -0.5022470
#> 6      1 -0.9151832 -0.15757508  0.1973134 -1.31430392  0.24629283  0.9822986
#>            7          8           9           10         11         12
#> 1  0.3711177 -0.7330453 -0.63540128 -0.505026336  1.1774047  0.4296625
#> 2 -0.8354430 -0.2938705  0.03111661  0.085975510  1.1477487  0.4023173
#> 3 -1.6674136 -1.0990787 -1.00007998 -0.947069513  0.3588851 -0.8873257
#> 4 -0.9068406 -0.8932880 -0.58649863 -0.683863298 -0.9053647 -0.4603737
#> 5 -0.6486205 -0.9952912 -0.44048783 -0.006790727 -0.5014230  0.1964397
#> 6 -0.5239507 -0.6489432 -0.53016266 -0.683924188  1.2357249 -0.2040576
#>           13         14         15         16         17          18         19
#> 1  0.4330098  0.2688759 -0.2456852  1.1310980  0.5107434 -1.49171524  0.5461195
#> 2 -0.4543559 -0.9603846 -1.4064577 -0.4327531 -0.6527845 -2.36113493  2.0883711
#> 3 -1.4776433 -1.7967629 -1.9298573 -0.7471992 -1.4625251 -1.10579922  1.4295638
#> 4 -0.6764329 -0.9565217 -0.3102249 -0.1738450 -0.6722008 -0.01761965 -0.3733054
#> 5 -0.8742965 -0.6689574 -0.3878315  0.1411472 -0.6709104 -1.22389933 -0.7333397
#> 6 -0.4211337 -0.1178140 -1.2466737  0.1826659 -0.4473040 -2.16117290 -0.1881705
#>            20      New_1       New_2      New_3       New_4      New_5
#> 1 -0.44740467 -1.1643327 -0.06426498 -1.0297345 -0.77352947 -0.5462547
#> 2  0.79348781 -2.5153602 -1.36237046 -2.0750023 -0.01391596 -1.2526368
#> 3  0.01670341 -1.5287918 -2.41782779 -1.2808399 -1.23474600 -2.2184071
#> 4 -0.31026289 -0.3648568 -0.20992439  0.1999604 -0.81191510 -0.8132315
#> 5 -0.38518119 -0.9982079 -0.59390388 -0.8748418 -0.61917135 -0.1031998
#> 6 -0.48971901 -1.9617002 -1.33014481 -1.4159325 -0.20824058 -1.3501761
```

We will hold out the new locations for parameter estimation.

``` r

data_old <- data_all[, 1:21]
data_new <- data_all[, -c(1:21)]
```

### Fitting Covariance Models

We will first fit pure spatial and temporal models, then the fully
symmetric model, and finally the general stationary model. First, we
will create an `mcgf_rs` object and calculate auto- and
cross-correlations.

``` r

data_mcgf <- mcgf_rs(data_old[, -1],
    locations = locations, longlat = TRUE,
    label = data_old[, 1]
)
#> `time` is not provided, assuming rows are equally spaced temporally.
data_mcgf <- add_acfs(data_mcgf, lag_max = lag)
data_mcgf <- add_ccfs(data_mcgf, lag_max = lag)

# If multiple cores are available
# data_mcgf <- add_ccfs(data_mcgf, lag_max = lag, ncores = 8)
```

Here `acfs` actually refers to the mean auto-correlations across the
stations for each time lag. To view the calculated `acfs`, we can run:

``` r

acfs(data_mcgf)
#> $acfs
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6075506 0.4076272 0.3805953 0.3702800 0.3586723 
#> 
#> $acfs_rs
#> $acfs_rs$`Regime 1`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.5953378 0.3991894 0.3814418 0.3714696 0.3639292 
#> 
#> $acfs_rs$`Regime 2`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6149751 0.4125410 0.3797082 0.3691936 0.3550066
```

Similarly, we can view the `ccfs` by the `ccfs` function. Here we will
present lag 1 regime-switching ccfs for the first 6 locations for the
two regimes.

``` r

# Regime 1
ccfs(data_mcgf)$ccfs_rs$`Regime 1`[1:6, 1:6, 2]
#>           1         2          3         4         5          6
#> 1 0.6047946 0.1236369 0.05451172 0.3713023 0.1363766 0.07120106
#> 2 0.3027375 0.6174314 0.25427725 0.4220703 0.4477891 0.23134358
#> 3 0.2293202 0.3347340 0.58222643 0.3837498 0.2952973 0.46537082
#> 4 0.5256692 0.1564803 0.12938351 0.6003370 0.1701296 0.11145141
#> 5 0.2602507 0.5444454 0.32282783 0.3751226 0.6132754 0.30494395
#> 6 0.1959693 0.3710529 0.53421219 0.3193299 0.3348082 0.56379838

# Regime 2
ccfs(data_mcgf)$ccfs_rs$`Regime 2`[1:6, 1:6, 2]
#>           1         2         3         4         5         6
#> 1 0.6172964 0.1575309 0.1723283 0.3893590 0.1096808 0.1688383
#> 2 0.2620954 0.6289374 0.3280881 0.3878912 0.4234188 0.3004865
#> 3 0.2811509 0.3158273 0.6163612 0.4355238 0.2379046 0.4932367
#> 4 0.5759686 0.1825519 0.2120796 0.6293679 0.1514506 0.1733759
#> 5 0.1509059 0.5799483 0.4202242 0.2930347 0.6249524 0.4089683
#> 6 0.2351466 0.3713949 0.5875119 0.3502719 0.3122138 0.6069874
```

#### Pure Spatial Model

The pure spatial model can be fitted using the `fit_base` function. The
results are actually obtained from the optimization function `nlminb`.
Since the base model is the same for the two regimes, we set
`rs = FALSE` to indicate this.

``` r

fit_spatial <- fit_base(
    data_mcgf,
    model_ls = "spatial",
    lag_ls = lag,
    par_init_ls = list(c(c = 0.000001, gamma = 0.5)),
    par_fixed_ls = list(c(nugget = 0)),
    rs = FALSE
)
fit_spatial[[1]]$fit
#> $par
#>          c      gamma 
#> 0.00390725 0.50000000 
#> 
#> $objective
#> [1] 13.20815
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 2
#> 
#> $evaluations
#> function gradient 
#>       11        7 
#> 
#> $message
#> [1] "both X-convergence and relative convergence (5)"
```

#### Pure Temporal Model

The pure temporal model can also be fitted by `fit_base`:

``` r

fit_temporal <- fit_base(
    data_mcgf,
    model = "temporal",
    lag_ls = lag,
    par_init_ls = list(c(a = 1, alpha = 0.5)),
    rs = FALSE
)
fit_temporal[[1]]$fit
#> $par
#>         a     alpha 
#> 0.7137588 0.3484975 
#> 
#> $objective
#> [1] 0.0202941
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 9
#> 
#> $evaluations
#> function gradient 
#>       11       26 
#> 
#> $message
#> [1] "relative convergence (4)"
```

We can store the fitted spatial and temporal models to `data_mcgf` using
`add_base`:

``` r

data_mcgf <- add_base(data_mcgf,
    fit_s = fit_spatial,
    fit_t = fit_temporal,
    sep = T
)
```

#### Separable Model

We can also calculate the WLS estimates all at once by fitting the
separable model:

``` r

fit_sep <- fit_base(
    data_mcgf,
    model_ls = "sep",
    lag_ls = lag,
    par_init_ls = list(c(c = 0.000001, gamma = 0.5, a = 1, alpha = 0.5)),
    par_fixed_ls = list(c(nugget = 0)),
    control = list(list(iter.max = 10000, eval.max = 10000)),
    rs = FALSE
)
fit_sep[[1]]$fit
#> $par
#>          c      gamma          a      alpha 
#> 0.00390725 0.50000000 0.72365072 0.39736171 
#> 
#> $objective
#> [1] 29.16779
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 11
#> 
#> $evaluations
#> function gradient 
#>       27       64 
#> 
#> $message
#> [1] "relative convergence (4)"
```

store the fully symmetric model as the base model and print the base
model:

``` r

data_mcgf <- add_base(data_mcgf, fit_base = fit_sep, old = TRUE)
model(data_mcgf, model = "base", old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Regime switching: FALSE
#> - Base model: sep 
#> - Parameters:
#>          c      gamma     nugget          a      alpha 
#> 0.00390725 0.50000000 0.00000000 0.71375884 0.34849745 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Regime switching: FALSE
#> - Base model: sep 
#> - Parameters:
#>          c      gamma          a      alpha     nugget 
#> 0.00390725 0.50000000 0.72365072 0.39736171 0.00000000 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb
```

The `old = TRUE` in `add_base` keeps the fitted pure spatial and
temporal models for records, and they are not used for any further
steps. It is recommended to keep the old model not only for
reproducibility, but to keep a history of fitted models.

#### Lagrangian Model

We will fit Lagrangian correlation functions to model the
regime-dependent prevailing winds:

``` r

fit_lagr_rs <- fit_lagr(data_mcgf,
    model_ls = list("lagr_tri"),
    par_init_ls = list(list(lambda = 0.1, v1 = 100, v2 = 100, k = 2))
)
lapply(fit_lagr_rs[1:2], function(x) x$fit)
#> $`Regime 1`
#> $`Regime 1`$par
#>     lambda         v1         v2          k 
#>  0.1702415 74.3744096 65.6736958  1.3422782 
#> 
#> $`Regime 1`$objective
#> [1] 12.76536
#> 
#> $`Regime 1`$convergence
#> [1] 0
#> 
#> $`Regime 1`$iterations
#> [1] 87
#> 
#> $`Regime 1`$evaluations
#> function gradient 
#>      116      513 
#> 
#> $`Regime 1`$message
#> [1] "relative convergence (4)"
#> 
#> 
#> $`Regime 2`
#> $`Regime 2`$par
#>     lambda         v1         v2          k 
#>  0.2167212 43.9009790 78.2938279  1.7852123 
#> 
#> $`Regime 2`$objective
#> [1] 15.68579
#> 
#> $`Regime 2`$convergence
#> [1] 0
#> 
#> $`Regime 2`$iterations
#> [1] 27
#> 
#> $`Regime 2`$evaluations
#> function gradient 
#>       35      128 
#> 
#> $`Regime 2`$message
#> [1] "relative convergence (4)"
```

Finally we will store this model to `data_mcgf` using `add_lagr` and
then print the final model:

``` r

data_mcgf <- add_lagr(data_mcgf, fit_lagr = fit_lagr_rs)
model(data_mcgf, old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 5, 5 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Regime switching: FALSE
#> - Base model: sep 
#> - Parameters:
#>          c      gamma     nugget          a      alpha 
#> 0.00390725 0.50000000 0.00000000 0.71375884 0.34849745 
#> 
#> - Fixed parameters:
#> nugget 
#>      0 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Regime switching: FALSE
#> - Base model: sep 
#> - Parameters:
#>          c      gamma          a      alpha     nugget 
#> 0.00390725 0.50000000 0.72365072 0.39736171 0.00000000 
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
#> - Regime switching: TRUE 
#> --------------------
#>       Regime 1 
#> --------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>     lambda         v1         v2          k 
#>  0.1702415 74.3744096 65.6736958  1.3422782 
#> 
#> - Fixed parameters:
#> NULL
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> --------------------
#>       Regime 2 
#> --------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>     lambda         v1         v2          k 
#>  0.2167212 43.9009790 78.2938279  1.7852123 
#> 
#> - Fixed parameters:
#> NULL
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb
```

### Simple Kriging for new locations

We provide functions for computing simple Kriging forecasts for new
locations. The associated function is `krige_new`, and users can either
supply the coordinates for the new locations or the distance matrices
for all locations.

``` r

krige_base_new <- krige_new(
    x = data_mcgf,
    locations_new = locations_new,
    model = "base",
    interval = TRUE
)
krige_stat_new <- krige_new(
    x = data_mcgf,
    locations_new = locations_new,
    model = "all",
    interval = TRUE
)
```

Next, we can compute the root mean square error (RMSE), mean absolute
error (MAE), and the realized percentage of observations falling outside
the 95% PI (POPI) for these models for the new locations.

#### RMSE

``` r

# RMSE
rmse_base <-
    sqrt(colMeans((data_new - krige_base_new$fit[, -c(1:20)])^2, na.rm = T))
rmse_stat <-
    sqrt(colMeans((data_new - krige_stat_new$fit[, -c(1:20)])^2, na.rm = T))
rmse <- c(
    "Base" = mean(rmse_base),
    "STAT" = mean(rmse_stat)
)
rmse
#>      Base      STAT 
#> 0.8117693 0.7622198
```

#### MAE

``` r

mae_base <- colMeans(abs(data_new - krige_base_new$fit[, -c(1:20)]), na.rm = T)
mae_stat <- colMeans(abs(data_new - krige_stat_new$fit[, -c(1:20)]), na.rm = T)
mae <- c(
    "Base" = mean(mae_base),
    "STAT" = mean(mae_stat)
)
mae
#>      Base      STAT 
#> 0.6451738 0.6061245
```

#### POPI

``` r

# POPI
popi_base <- colMeans(
    data_new < krige_base_new$lower[, -c(1:20)] | data_new > krige_base_new$upper[, -c(1:20)],
    na.rm = T
)
popi_stat <- colMeans(
    data_new < krige_stat_new$lower[, -c(1:20)] | data_new > krige_stat_new$upper[, -c(1:20)],
    na.rm = T
)
popi <- c(
    "Base" = mean(popi_base),
    "STAT" = mean(popi_stat)
)
popi
#>       Base       STAT 
#> 0.04560561 0.03811812
```

## References

Gneiting, Tilmann. 2002. “Nonseparable, Stationary Covariance Functions
for Space–Time Data.” *Journal of the American Statistical Association*
97 (458): 590–600. <https://doi.org/10.1198/016214502760047113>.

Jia, Tianxia, and Deniz Sezer. 2025. “Improving Short-Term Wind Speed
Forecasts Using Regime-Switching Spatio-Temporal Covariance Models.”
*Journal of Renewable and Sustainable Energy* 17 (5): 053304.
<https://doi.org/10.1063/5.0285012>.
