# Calculate standard deviation for each location under each regime.

Calculate standard deviation for each location under each regime.

## Usage

``` r
sd_rs(x, label)
```

## Arguments

- x:

  A `data.frame` or `matrix`.

- label:

  A vector of regime labels. Its length must be the same as the number
  of rows in `x`.

## Value

A list of standard deviations for each regime.

## Examples

``` r
set.seed(123)
x <- matrix(rnorm(200), nrow = 100)
label <- sample(1:2, 100, replace = TRUE)
sd_rs(x, label = factor(label))
#> $`Regime 1`
#> [1] 1.0091067 0.9136253
#> 
#> $`Regime 2`
#> [1] 0.7816577 1.0241517
#> 
```
