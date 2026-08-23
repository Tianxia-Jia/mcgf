# Extract, calculate, or assign mean auto-correlations for an `mcgf` or `mcgf_rs` object

Extract, calculate, or assign mean auto-correlations for an `mcgf` or
`mcgf_rs` object

## Usage

``` r
# S3 method for class 'mcgf'
acfs(x, lag_max, replace = FALSE, ...)

# S3 method for class 'mcgf_rs'
acfs(x, lag_max, replace = FALSE, ...)

acfs(x) <- value

add_acfs(x, lag_max, ...)
```

## Arguments

- x:

  An `mcgf` or `mcgf_rs` object.

- lag_max:

  Maximum lag at which to calculate the acf.

- replace:

  Logical; if TRUE, `acfs` are recalculated.

- ...:

  Additional parameters or attributes.

- value:

  A Vector of mean of auto-correlations for time lags starting from 0.

## Value

[`acfs()`](https://tianxia-jia.github.io/mcgf/reference/acfs.md) returns
(regime-switching) mean auto-correlations. `add_acfs()` returns the same
object with additional attributes of (regime-switching) mean
auto-correlations.

## Details

For `mcgf` objects,
[`acfs()`](https://tianxia-jia.github.io/mcgf/reference/acfs.md)
computes mean auto-correlations for each time lag across locations. The
output is a vector of acfs.

For `mcgf_rs` objects,
[`acfs()`](https://tianxia-jia.github.io/mcgf/reference/acfs.md)
computes regime-switching mean auto-correlations for each time lag
across locations. The output is a list of vectors of acfs, where each
vector in the list corresponds to the acfs for a regime.

`acfs<-` assigns `acfs` to `x`.

`add_acfs()` adds `acfs` to `x`.

## Examples

``` r
# Calculate acfs for 'sim1'
data(sim1)
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
acfs(sim1_mcgf, lag_max = 5)
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6109553 0.3140734 0.2401463 0.1781788 0.1351682 

# Add acfs to 'sim1_mcgf'
sim1_mcgf <- add_acfs(sim1_mcgf, lag_max = 5)
print(sim1_mcgf, "acfs")
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6109553 0.3140734 0.2401463 0.1781788 0.1351682 

# Calculate acfs for 'sim2'
data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
acfs(sim2_mcgf, lag_max = 5)
#> $acfs
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.7195535 0.5364966 0.4539327 0.3834479 0.3308474 
#> 
#> $acfs_rs
#> $acfs_rs$`Regime 1`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6903615 0.6021888 0.5731002 0.4950374 0.4379770 
#> 
#> $acfs_rs$`Regime 2`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.7558502 0.4491320 0.2978408 0.2379951 0.1911654 
#> 
#> 

# Add acfs to 'sim2_mcgf'
sim2_mcgf <- add_acfs(sim2_mcgf, lag_max = 5)
print(sim2_mcgf, "acfs")
#> $acfs
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.7195535 0.5364966 0.4539327 0.3834479 0.3308474 
#> 
#> $acfs_rs
#> $acfs_rs$`Regime 1`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.6903615 0.6021888 0.5731002 0.4950374 0.4379770 
#> 
#> $acfs_rs$`Regime 2`
#>      lag0      lag1      lag2      lag3      lag4      lag5 
#> 1.0000000 0.7558502 0.4491320 0.2978408 0.2379951 0.1911654 
#> 
#> 
```
