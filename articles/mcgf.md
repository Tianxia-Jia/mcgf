# Steps for fitting an \`mcgf\` object

## Introduction

The `mcgf` package contains useful functions to simulate Markov chain
Gaussian fields (MCGF) and regime-switching Markov chain Gaussian fields
(RS-MCGF) with covariance structures of the Gneiting class (Gneiting
2002). It also provides useful tools to estimate the parameters by
weighted least squares (WLS) and maximum likelihood estimation (MLE).
The `mcgf` function can be used to fit covariance models and obtain
Kriging forecasts. A typical workflow for fitting a non-regime-switching
mcgf is given below.

1.  Create an `mcgf` object by providing a dataset and the corresponding
    locations/distances
2.  Calculate auto- and cross-correlations.
3.  Fit the base covariance model which is a fully symmetric model.
4.  Fit the Lagrangian model to account for asymmetry, if necessary.
5.  Obtain Kriging forecasts.
6.  Obtain Kriging forecasts for new locations given their coordinates.

We will demonstrate the use of `mcgf` by an example below.

## Irish Wind

The Ireland wind data contains daily wind speeds for 1961-1978 at 11
synoptic meteorological stations in the Republic of Ireland. It is
available in the `gstat` package and is also imported to `mcgf`. To view
the details on this dataset, run
[`help(wind)`](https://tianxia-jia.github.io/mcgf/reference/wind.md).
The object `wind` contains the wind speeds data as well as the locations
of the stations.

``` r

library(mcgf)
data(wind)
head(wind$data)
#>         time      VAL      BEL       CLA      SHA      ROC      BIR      MUL
#> 1 1961-01-01 7.696089 9.517222 5.2730556 7.181644 7.737244 5.077567 5.571433
#> 2 1961-01-02 8.683822 9.023356 5.1650222 6.492289 7.567478 3.945789 5.036411
#> 3 1961-01-03 8.683822 6.559167 4.1361333 5.746344 9.517222 3.174122 4.372778
#> 4 1961-01-04 3.410767 2.808867 0.9208556 2.335578 5.442822 1.481600 2.999211
#> 5 1961-01-05 6.816389 6.646622 3.3644667 5.509700 6.857544 4.223589 5.617733
#> 6 1961-01-06 4.177289 4.177289 2.2738444 2.762567 6.795811 2.315000 3.688567
#>        MAL      KIL      CLO      DUB
#> 1 7.737244 4.779189 6.471711 7.032456
#> 2 7.114767 3.343889 4.974678 5.916111
#> 3 6.538589 5.211322 3.945789 5.787500
#> 4 5.597156 2.356156 3.024933 4.439656
#> 5 6.085878 3.174122 5.319356 6.132178
#> 6 6.775233 3.431344 3.858333 5.489122
wind$locations
#>            lon      lat
#> VAL -10.250000 51.93333
#> BEL -10.000000 54.23333
#> CLA  -8.983333 53.71667
#> SHA  -8.916667 52.70000
#> ROC  -8.250000 51.80000
#> BIR  -7.883333 53.08333
#> MUL  -7.366667 53.53333
#> MAL  -7.333333 55.36667
#> KIL  -7.266667 52.66667
#> CLO  -7.233333 54.18333
#> DUB  -6.250000 53.43333
```

### Data De-trending

To fit covariance models on `wind`, the first step is to load the data
set and de-trend the data to have a mean-zero time series. We will
follow the data de-trending procedure in Gneiting et al. (2006).

1.  We start with removing the leap day and take the square-root
    transformation. To do this we need to use the `lubridate` package.

``` r

# install.packages("lubridate")
library(lubridate)
ind_leap <- month(wind$data$time) == 2 & day(wind$data$time) == 29
wind_de <- wind$data[!ind_leap, ]
wind_de[, -1] <- sqrt(wind_de[, -1])
```

2.  Next, we split the data into training and test datasets, where
    training contains data in 1961-1970 and test spans 1971-1978.

``` r

is_train <- year(wind_de$time) <= 1970
wind_train <- wind_de[is_train, ]
wind_test <- wind_de[!is_train, ]
```

3.  We will estimate and remove the annual trend for the training
    dataset. The annual trend is the daily averages across the training
    years. We will use the `dplyr` package for this task.

``` r

# install.packages("dplyr")
library(dplyr)
# Estimate annual trend
avg_train <- tibble(
    month = month(wind_train$time),
    day = day(wind_train$time),
    ws = rowMeans(wind_train[, -1])
)
trend_train <- avg_train %>%
    group_by(month, day) %>%
    summarise(trend = mean(ws))

# Subtract annual trend
trend_train2 <- left_join(avg_train, trend_train, by = c("month", "day"))$trend
wind_train[, -1] <- wind_train[, -1] - trend_train2
```

4.  To obtain a mean-zero time series, we will subtract the station-wise
    mean for each station.

``` r

# Subtract station-wise mean
mean_train <- colMeans(wind_train[, -1])
wind_train[, -1] <- sweep(wind_train[, -1], 2, mean_train)
wind_trend <- list(
    annual = as.data.frame(trend_train),
    mean = mean_train
)
```

5.  Finally, we will subtract the annual trend and station-wise mean
    from the test dataset.

``` r

avg_test <- tibble(
    month = month(wind_test$time),
    day = day(wind_test$time)
)
trend_train3 <-
    left_join(avg_test, trend_train, by = c("month", "day"))$trend
wind_test[, -1] <- wind_test[, -1] - trend_train3
wind_test[, -1] <- sweep(wind_test[, -1], 2, mean_train)
```

### Fitting Covariance Models

We will first fit pure spatial and temporal models, then the separable
model, and finally the general stationary model. First, we will create
an `mcgf` object and calculate auto- and cross-correlations.

``` r

wind_mcgf <- mcgf(wind_train[, -1], locations = wind$locations, longlat = TRUE)
#> `time` is not provided, assuming rows are equally spaced temporally.
wind_mcgf <- add_acfs(wind_mcgf, lag_max = 3)
wind_mcgf <- add_ccfs(wind_mcgf, lag_max = 3)
```

Here `acfs` actually refers to the mean auto-correlations across the
stations for each time lag. To view the calculated `acfs`, we can run:

``` r

acfs(wind_mcgf)
#>      lag0      lag1      lag2      lag3 
#> 1.0000000 0.5072739 0.2382621 0.1643998
```

Similarly, we can view the `ccfs` by:

``` r

ccfs(wind_mcgf)
#> , , lag0
#> 
#>           VAL       BEL       CLA       SHA       ROC       BIR       MUL
#> VAL 1.0000000 0.7076428 0.7588797 0.8299603 0.7962082 0.7719814 0.6688609
#> BEL 0.7076428 1.0000000 0.8471918 0.7432194 0.5795413 0.7650007 0.7371354
#> CLA 0.7588797 0.8471918 1.0000000 0.8520358 0.6994680 0.8777462 0.8410719
#> SHA 0.8299603 0.7432194 0.8520358 1.0000000 0.7960579 0.8865147 0.8239009
#> ROC 0.7962082 0.5795413 0.6994680 0.7960579 1.0000000 0.7592625 0.7305102
#> BIR 0.7719814 0.7650007 0.8777462 0.8865147 0.7592625 1.0000000 0.8841425
#> MUL 0.6688609 0.7371354 0.8410719 0.8239009 0.7305102 0.8841425 1.0000000
#> MAL 0.5153690 0.7087588 0.6818106 0.6121929 0.5357199 0.6411302 0.7375539
#> KIL 0.7333467 0.6487696 0.7772730 0.8236153 0.8338255 0.8406781 0.8414494
#> CLO 0.6718720 0.7697673 0.8511102 0.7759129 0.6876764 0.8436560 0.8699373
#> DUB 0.5967529 0.6498491 0.7420593 0.7361857 0.6862624 0.7812147 0.8609285
#>           MAL       KIL       CLO       DUB
#> VAL 0.5153690 0.7333467 0.6718720 0.5967529
#> BEL 0.7087588 0.6487696 0.7697673 0.6498491
#> CLA 0.6818106 0.7772730 0.8511102 0.7420593
#> SHA 0.6121929 0.8236153 0.7759129 0.7361857
#> ROC 0.5357199 0.8338255 0.6876764 0.6862624
#> BIR 0.6411302 0.8406781 0.8436560 0.7812147
#> MUL 0.7375539 0.8414494 0.8699373 0.8609285
#> MAL 1.0000000 0.5892218 0.7568132 0.7145366
#> KIL 0.5892218 1.0000000 0.7934818 0.7770512
#> CLO 0.7568132 0.7934818 1.0000000 0.8161496
#> DUB 0.7145366 0.7770512 0.8161496 1.0000000
#> 
#> , , lag1
#> 
#>           VAL       BEL       CLA       SHA       ROC       BIR       MUL
#> VAL 0.4739681 0.4078097 0.3734029 0.4028977 0.3406302 0.3871317 0.3171416
#> BEL 0.3843919 0.5310138 0.4208376 0.3899411 0.2819817 0.4000168 0.3695143
#> CLA 0.4164508 0.4814516 0.4993218 0.4450647 0.3244728 0.4521809 0.4020938
#> SHA 0.4719522 0.4646048 0.4689832 0.5137855 0.3883962 0.4677375 0.4100365
#> ROC 0.4706594 0.3793417 0.3852120 0.4285917 0.4442973 0.4113073 0.3590478
#> BIR 0.4725125 0.4964189 0.4986726 0.4927316 0.3919890 0.5409202 0.4557624
#> MUL 0.4438395 0.5082231 0.5039360 0.4898975 0.4056061 0.5097809 0.5228395
#> MAL 0.3550724 0.4616666 0.4118390 0.3918386 0.3360198 0.3792364 0.4182625
#> KIL 0.4867586 0.4607927 0.4866086 0.4964729 0.4299541 0.4930438 0.4430179
#> CLO 0.4390459 0.5117641 0.5057839 0.4686427 0.3853159 0.4904330 0.4549374
#> DUB 0.4333070 0.4698080 0.4829445 0.4819004 0.4349093 0.4860929 0.4858079
#>           MAL       KIL       CLO       DUB
#> VAL 0.2712570 0.2917876 0.3248601 0.2900024
#> BEL 0.3554572 0.3012227 0.3806926 0.3179568
#> CLA 0.3289562 0.3458138 0.4112028 0.3569920
#> SHA 0.3203410 0.3834180 0.4088611 0.3736715
#> ROC 0.2950458 0.3540329 0.3479210 0.3484134
#> BIR 0.3408605 0.3973215 0.4505508 0.3984836
#> MUL 0.4151149 0.4232129 0.4664284 0.4434980
#> MAL 0.5296408 0.3105146 0.4036000 0.4002747
#> KIL 0.3227371 0.4642266 0.4387080 0.3981810
#> CLO 0.4032414 0.4076724 0.5190745 0.4289406
#> DUB 0.4428008 0.4300411 0.4765407 0.5409250
#> 
#> , , lag2
#> 
#>           VAL       BEL       CLA       SHA       ROC       BIR       MUL
#> VAL 0.2011319 0.1973159 0.1630851 0.1871285 0.1489858 0.1943355 0.1492425
#> BEL 0.1713197 0.2588323 0.1891480 0.1936055 0.1240906 0.2020276 0.1806373
#> CLA 0.1585236 0.2132206 0.2247423 0.2000905 0.1184259 0.2115274 0.1763905
#> SHA 0.1938118 0.2298367 0.2157028 0.2391773 0.1589136 0.2265611 0.1844087
#> ROC 0.1765714 0.1716184 0.1421018 0.1764547 0.1696931 0.1742268 0.1451217
#> BIR 0.2020194 0.2401048 0.2278565 0.2262318 0.1628422 0.2765692 0.2136481
#> MUL 0.1862142 0.2487854 0.2182151 0.2159567 0.1694341 0.2448788 0.2510983
#> MAL 0.1555252 0.2110852 0.1623937 0.1691266 0.1462008 0.1585537 0.1853651
#> KIL 0.1743715 0.2068322 0.1970328 0.2091868 0.1495429 0.2192173 0.1882474
#> CLO 0.1835787 0.2349159 0.2190671 0.2138391 0.1511675 0.2355708 0.2042911
#> DUB 0.1856084 0.2176080 0.2022041 0.2148639 0.1950841 0.2234817 0.2150465
#>           MAL       KIL       CLO       DUB
#> VAL 0.1427551 0.1382869 0.1666310 0.1535865
#> BEL 0.1754907 0.1388237 0.1941650 0.1629718
#> CLA 0.1419690 0.1537041 0.1979937 0.1674888
#> SHA 0.1503395 0.1803397 0.2022438 0.1764521
#> ROC 0.1458883 0.1373046 0.1488279 0.1654739
#> BIR 0.1438303 0.1844782 0.2197392 0.1920556
#> MUL 0.1950942 0.1905402 0.2183714 0.2084112
#> MAL 0.2636079 0.1107242 0.1729375 0.1921164
#> KIL 0.1334084 0.1931797 0.1982144 0.1790683
#> CLO 0.1715671 0.1793731 0.2624348 0.2135477
#> DUB 0.2180694 0.1899728 0.2270066 0.2804169
#> 
#> , , lag3
#> 
#>            VAL       BEL        CLA        SHA        ROC        BIR        MUL
#> VAL 0.15186704 0.1424697 0.12263507 0.12486035 0.11217575 0.14386635 0.09991865
#> BEL 0.12223738 0.1854410 0.13490665 0.13169403 0.08728794 0.14207057 0.12397596
#> CLA 0.10549441 0.1408343 0.16084299 0.13066846 0.07914442 0.13859706 0.10865229
#> SHA 0.12004043 0.1511578 0.14232764 0.14961877 0.09888111 0.14670946 0.11247087
#> ROC 0.11051322 0.1022188 0.08007433 0.09282562 0.10856180 0.10063755 0.07957425
#> BIR 0.13309457 0.1594322 0.15267585 0.13773090 0.10595023 0.19194636 0.13771463
#> MUL 0.11364183 0.1682815 0.14195183 0.12747872 0.10407473 0.15372629 0.16630054
#> MAL 0.09065312 0.1431289 0.10244646 0.10197080 0.09770836 0.08625211 0.11386343
#> KIL 0.09511038 0.1179299 0.11535092 0.11479319 0.07924395 0.12369722 0.10034132
#> CLO 0.12004487 0.1584521 0.15581710 0.13825387 0.10650927 0.15601087 0.13152273
#> DUB 0.10526163 0.1331954 0.11965370 0.11622550 0.11602725 0.12380225 0.11590933
#>            MAL        KIL        CLO       DUB
#> VAL 0.11038186 0.09840820 0.12881956 0.1131093
#> BEL 0.11895940 0.09539470 0.13930407 0.1065024
#> CLA 0.08566429 0.10250903 0.13807158 0.1059626
#> SHA 0.10015324 0.11889068 0.14364630 0.1110210
#> ROC 0.10700118 0.08017778 0.09975703 0.1043980
#> BIR 0.09283539 0.12122167 0.15906405 0.1304590
#> MUL 0.13033344 0.11347989 0.14872192 0.1286318
#> MAL 0.19778972 0.05899522 0.10990283 0.1253038
#> KIL 0.07627083 0.11333371 0.12609113 0.1024267
#> CLO 0.10926382 0.12014806 0.20136056 0.1459052
#> DUB 0.14141118 0.10085759 0.14688397 0.1813353
```

#### Pure Spatial Model

The pure spatial model can be fitted using the `fit_base` function. The
results are actually obtained from the optimization function `nlminb`.

``` r

fit_spatial <- fit_base(
    wind_mcgf,
    model = "spatial",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001),
    par_fixed = c(gamma = 0.5)
)
fit_spatial$fit
#> $par
#>           c      nugget 
#> 0.001326192 0.049136042 
#> 
#> $objective
#> [1] 1.742162
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 9
#> 
#> $evaluations
#> function gradient 
#>       23       24 
#> 
#> $message
#> [1] "relative convergence (4)"
```

Here we set `gamma` to be 0.5 and it is not estimated along with `c` or
`nugget`. By default `mcgf` provides two optimization functions:
`nlminb` and `optim`. Other optimization functions are also supported as
long as their first two arguments are initial values for the parameters
and a function to be minimized respectively (same as that of `optim` and
`nlminb`). Also, if the argument names for upper and lower bounds are
not `upper` or `lower`, we can create a simple wrapper to “change” them.

``` r

library(Rsolnp)
solnp2 <- function(pars, fun, lower, upper, ...) {
    solnp(pars, fun, LB = lower, UB = upper, ...)
}
fit_spatial2 <- fit_base(
    wind_mcgf,
    model = "spatial",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001),
    par_fixed = c(gamma = 0.5),
    optim_fn = "other",
    other_optim_fn = "solnp2"
)
#> 
#> Iter: 1 fn: 1.7422    Pars:  0.001326 0.049164
#> Iter: 2 fn: 1.7422    Pars:  0.001326 0.049164
#> solnp--> Completed in 2 iterations
fit_spatial2$fit
#> $pars
#>           c      nugget 
#> 0.001325958 0.049163635 
#> 
#> $convergence
#> [1] 0
#> 
#> $values
#> [1] 2.539010 1.742162 1.742162
#> 
#> $lagrange
#> [1] 0
#> 
#> $hessian
#>            [,1]       [,2]
#> [1,] 58504276.4 391217.107
#> [2,]   391217.1   3319.562
#> 
#> $ineqx0
#> NULL
#> 
#> $nfuneval
#> [1] 79
#> 
#> $outer.iter
#> [1] 2
#> 
#> $elapsed
#> Time difference of 0.008444071 secs
#> 
#> $vscale
#> [1] 1 1 1
```

#### Pure Temporal Model

The pure temporal model can also be fitted by `fit_base`:

``` r

fit_temporal <- fit_base(
    wind_mcgf,
    model = "temporal",
    lag = 3,
    par_init = c(a = 0.5, alpha = 0.5)
)
fit_temporal$fit
#> $par
#>         a     alpha 
#> 0.9774157 0.8053363 
#> 
#> $objective
#> [1] 0.0006466064
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 12
#> 
#> $evaluations
#> function gradient 
#>       15       31 
#> 
#> $message
#> [1] "both X-convergence and relative convergence (5)"
```

Before fitting the fully symmetric model, we need to store the fitted
spatial and temporal models to `wind_mcgf` using `add_base`:

``` r

wind_mcgf <- add_base(wind_mcgf,
    fit_s = fit_spatial,
    fit_t = fit_temporal,
    sep = T
)
```

#### Separable Model

We can also fit the pure spatial and temporal models all at once by
fitting a separable model:

``` r

fit_sep <- fit_base(
    wind_mcgf,
    model = "sep",
    lag = 3,
    par_init = c(nugget = 0.1, c = 0.001, a = 0.5, alpha = 0.5),
    par_fixed = c(gamma = 0.5)
)
fit_sep$fit
#> $par
#>           c      nugget           a       alpha 
#> 0.001309486 0.050562930 0.877171418 0.869117992 
#> 
#> $objective
#> [1] 2.838695
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 25
#> 
#> $evaluations
#> function gradient 
#>       54      121 
#> 
#> $message
#> [1] "relative convergence (4)"
```

One can also store this model to `wind_mcgf`, but to be consistent with
Gneiting et al. (2006), we will just use the estimates from the pure
spatial and temporal models to estimate the interaction parameter
`beta`.

#### Fully Symmetric Model

When holding other parameters constant, the parameter `beta` can be
estimated by:

``` r

par_sep <- c(fit_spatial$fit$par, fit_temporal$fit$par, gamma = 0.5)
fit_fs <-
    fit_base(
        wind_mcgf,
        model = "fs",
        lag = 3,
        par_init = c(beta = 0),
        par_fixed = par_sep
    )
fit_fs$fit
#> $par
#>      beta 
#> 0.6232458 
#> 
#> $objective
#> [1] 2.858384
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 5
#> 
#> $evaluations
#> function gradient 
#>        6        7 
#> 
#> $message
#> [1] "relative convergence (4)"
```

At this stage, we have fitted the base model, and we will store the
fully symmetric model as the base model and print the base model:

``` r

wind_mcgf <- add_base(wind_mcgf, fit_base = fit_fs, old = TRUE)
model(wind_mcgf, model = "base", old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 3 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Base-old model: sep 
#> - Parameters:
#>           c      nugget       gamma           a       alpha 
#> 0.001326192 0.049136042 0.500000000 0.977415717 0.805336266 
#> 
#> - Fixed parameters:
#> gamma 
#>   0.5 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Base model: fs 
#> - Parameters:
#>        beta           c      nugget           a       alpha       gamma 
#> 0.623245823 0.001326192 0.049136042 0.977415717 0.805336266 0.500000000 
#> 
#> - Fixed parameters:
#>           c      nugget           a       alpha       gamma 
#> 0.001326192 0.049136042 0.977415717 0.805336266 0.500000000 
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

We will fit a Lagrangian correlation function to model the westerly
wind:

``` r

fit_lagr <- fit_lagr(wind_mcgf,
    model = "lagr_tri",
    par_init = c(v1 = 200, lambda = 0.1),
    par_fixed = c(v2 = 0, k = 2)
)
fit_lagr$fit
#> $par
#>       lambda           v1 
#>   0.05595016 212.29485854 
#> 
#> $objective
#> [1] 2.662327
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 18
#> 
#> $evaluations
#> function gradient 
#>       21       54 
#> 
#> $message
#> [1] "relative convergence (4)"
```

Finally we will store this model to `wind_mcgf` using `add_lagr` and
then print the final model:

``` r

wind_mcgf <- add_lagr(wind_mcgf, fit_lagr = fit_lagr)
model(wind_mcgf, old = TRUE)
#> ----------------------------------------
#>                  Model
#> ----------------------------------------
#> - Time lag: 3 
#> - Scale of time lag: 1 
#> - Forecast horizon: 1 
#> ----------------------------------------
#>             Old - not in use
#> ----------------------------------------
#> - Base-old model: sep 
#> - Parameters:
#>           c      nugget       gamma           a       alpha 
#> 0.001326192 0.049136042 0.500000000 0.977415717 0.805336266 
#> 
#> - Fixed parameters:
#> gamma 
#>   0.5 
#> 
#> - Parameter estimation method: wls wls 
#> 
#> - Optimization function: nlminb nlminb 
#> ----------------------------------------
#>                  Base
#> ----------------------------------------
#> - Base model: fs 
#> - Parameters:
#>        beta           c      nugget           a       alpha       gamma 
#> 0.623245823 0.001326192 0.049136042 0.977415717 0.805336266 0.500000000 
#> 
#> - Fixed parameters:
#>           c      nugget           a       alpha       gamma 
#> 0.001326192 0.049136042 0.977415717 0.805336266 0.500000000 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb 
#> ----------------------------------------
#>               Lagrangian
#> ----------------------------------------
#> - Lagrangian model: lagr_tri 
#> - Parameters:
#>       lambda           v1           v2            k 
#>   0.05595016 212.29485854   0.00000000   2.00000000 
#> 
#> - Fixed parameters:
#> v2  k 
#>  0  2 
#> 
#> - Parameter estimation method: wls 
#> 
#> - Optimization function: nlminb
```

### Simple Kriging

For the test dataset, we can obtain the simple kriging forecasts and
intervals for the empirical model, base model, and general stationary
model using the `krige` function:

``` r

krige_emp <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "empirical",
    interval = TRUE
)
krige_base <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "base",
    interval = TRUE
)
krige_stat <- krige(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    model = "all",
    interval = TRUE
)
```

Next, we can compute the root mean square error (RMSE), mean absolute
error (MAE), and the realized percentage of observations falling outside
the 95% PI (POPI) for these models on the test dataset.

#### RMSE

``` r

# RMSE
rmse_emp <- sqrt(colMeans((wind_test[, -1] - krige_emp$fit)^2, na.rm = T))
rmse_base <- sqrt(colMeans((wind_test[, -1] - krige_base$fit)^2, na.rm = T))
rmse_stat <- sqrt(colMeans((wind_test[, -1] - krige_stat$fit)^2, na.rm = T))
rmse <- c(
    "Empirical" = mean(rmse_emp),
    "Base" = mean(rmse_base),
    "STAT" = mean(rmse_stat)
)
rmse
#> Empirical      Base      STAT 
#> 0.4778542 0.4890462 0.4852207
```

#### MAE

``` r

mae_emp <- colMeans(abs(wind_test[, -1] - krige_emp$fit), na.rm = T)
mae_base <- colMeans(abs(wind_test[, -1] - krige_base$fit), na.rm = T)
mae_stat <- colMeans(abs(wind_test[, -1] - krige_stat$fit), na.rm = T)
mae <- c(
    "Empirical" = mean(mae_emp),
    "Base" = mean(mae_base),
    "STAT" = mean(mae_stat)
)
mae
#> Empirical      Base      STAT 
#> 0.3776108 0.3890038 0.3851617
```

#### POPI

``` r

# POPI
popi_emp <- colMeans(
    wind_test[, -1] < krige_emp$lower | wind_test[, -1] > krige_emp$upper,
    na.rm = T
)
popi_base <- colMeans(
    wind_test[, -1] < krige_base$lower | wind_test[, -1] > krige_base$upper,
    na.rm = T
)
popi_stat <- colMeans(
    wind_test[, -1] < krige_stat$lower | wind_test[, -1] > krige_stat$upper,
    na.rm = T
)
popi <- c(
    "Empirical" = mean(popi_emp),
    "Base" = mean(popi_base),
    "STAT" = mean(popi_stat)
)
popi
#>  Empirical       Base       STAT 
#> 0.06962321 0.06473026 0.06326550
```

### Simple Kriging for new locations

We provide functions for computing simple Kriging forecasts for new
locations. The associated function is `krige_new`, and users can either
supply the coordinates for the new locations or the distance matrices
for all locations. Hypothetically, suppose a new location is at
($`9^\circ`$W, $`52^\circ`$N), then its kriging forecasts along with
forecasts for the rest of the stations for the general stationary model
can be obtained as follows.

``` r

krige_stat_new <- krige_new(
    x = wind_mcgf,
    newdata = wind_test[, -1],
    locations_new = c(-9, 52),
    model = "all",
    interval = TRUE
)
head(krige_stat_new$fit)
#>             VAL        BEL        CLA        SHA        ROC        BIR
#> 3653         NA         NA         NA         NA         NA         NA
#> 3654         NA         NA         NA         NA         NA         NA
#> 3655         NA         NA         NA         NA         NA         NA
#> 3656 -0.1447040 -0.3018214 -0.5758966 -0.4693314 -0.1046231 -0.5499500
#> 3657  0.5307650  0.1531176  0.2150980  0.3479844  0.4910715  0.3157901
#> 3658  0.1632995  0.1648983  0.1734723  0.1751080  0.1567910  0.3091076
#>             MUL        MAL        KIL        CLO        DUB      New_1
#> 3653         NA         NA         NA         NA         NA         NA
#> 3654         NA         NA         NA         NA         NA         NA
#> 3655         NA         NA         NA         NA         NA         NA
#> 3656 -0.5196208 -0.2799845 -0.4627520 -0.7397430 -0.4240240 -0.4152704
#> 3657  0.2576929  0.1404444  0.1013024  0.1157068  0.2114298  0.8410688
#> 3658  0.3274626  0.3504309  0.1426461  0.1006568  0.4521153  0.2792891
```

Below is a time-series plot for the first 100 forecasts for `New_1`:

``` r

plot.ts(krige_stat_new$fit[1:100, 12], ylab = "Wind Speed for New_1")
```

![](mcgf_files/figure-html/krige_plot-1.png)

## References

Gneiting, Tilmann. 2002. “Nonseparable, Stationary Covariance Functions
for Space–Time Data.” *Journal of the American Statistical Association*
97 (458): 590–600. <https://doi.org/10.1198/016214502760047113>.

Gneiting, Tilmann, Marc G. Genton, and Peter Guttorp. 2006.
*Geostatistical Space-Time Models, Stationarity, Separability and Full
Symmetry*. Department of Statistics, University of Washington;
<https://stat.uw.edu/sites/default/files/files/reports/2005/tr475.pdf>.
