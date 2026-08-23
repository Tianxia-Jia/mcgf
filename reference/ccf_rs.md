# Calculate regime-switching cross-correlation

Calculate regime-switching cross-correlation

## Usage

``` r
ccf_rs(x, y, label, lag_max)
```

## Arguments

- x, y:

  A univariate numeric time series.

- label:

  A factor of regime labels.

- lag_max:

  Maximum lag at which to calculate the ccf.

## Value

Cross-correlations for each group in `label`.

## Examples

``` r
set.seed(123)
x <- rnorm(100)
y <- rnorm(100)
label <- sample(1:2, 100, replace = TRUE)
ccf_rs(x, y, label = factor(label), lag_max = 3)
#> $`Regime 1`
#>          -3          -2          -1           0           1           2 
#>  0.07046131 -0.18742586 -0.07060330  0.01541658 -0.13212906  0.11329054 
#>           3 
#> -0.03214489 
#> 
#> $`Regime 2`
#>         -3         -2         -1          0          1          2          3 
#> -0.0191234 -0.3514831  0.1303284 -0.1254412  0.1901238  0.1122576  0.1892601 
#> 
```
