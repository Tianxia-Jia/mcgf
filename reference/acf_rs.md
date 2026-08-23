# Calculate regime-switching auto-correlation

Calculate regime-switching auto-correlation

## Usage

``` r
acf_rs(x, label, lag_max, demean = TRUE)
```

## Arguments

- x:

  A univariate numeric time series.

- label:

  A factor of regime labels.

- lag_max:

  Maximum lag at which to calculate the acf.

- demean:

  Logical. Should the covariances be about the sample means?

## Value

Mean auto-correlations for each group in `label`.

## Examples

``` r
set.seed(123)
x <- rnorm(100)
label <- sample(1:2, 100, replace = TRUE)
acf_rs(x, label = factor(label), lag_max = 3)
#> $`Regime 1`
#>        lag0        lag1        lag2        lag3 
#>  1.00000000 -0.06719448 -0.28622386  0.05784944 
#> 
#> $`Regime 2`
#>       lag0       lag1       lag2       lag3 
#> 1.00000000 0.02239601 0.08759688 0.25413972 
#> 
```
