# Convert correlation to covariance

Convert correlation to covariance

## Usage

``` r
cor2cov(V, sd, empirical = FALSE)

cor2cov_ar(V, sd, empirical = FALSE)
```

## Arguments

- V:

  A correlation matrix, usually positive semi-definite.

- sd:

  A vector of standard deviations.

- empirical:

  Logical; TRUE if V is empirical correlation.

## Value

A covariance matrix or covariance array.

## Details

`cor2cov` converts a matrix. `cor2cov_ar` converts an 3-D array.

## Examples

``` r
V <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
sd <- 1:2
cor2cov(V, sd)
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    4

V_ar <- array(c(1, 0.5, 0.5, 1), dim = c(2, 2, 2))
cor2cov_ar(V_ar, sd)
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    4
#> 
#> , , 2
#> 
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    4
#> 
```
