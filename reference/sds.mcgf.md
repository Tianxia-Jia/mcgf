# Extract, calculate, or assign standard deviations for an `mcgf` or `mcgf_rs` object.

Extract, calculate, or assign standard deviations for an `mcgf` or
`mcgf_rs` object.

## Usage

``` r
# S3 method for class 'mcgf'
sds(x, ...)

# S3 method for class 'mcgf_rs'
sds(x, replace = FALSE, ...)

sds(x) <- value
```

## Arguments

- x:

  An `mcgf` or `mcgf_rs` object.

- ...:

  Additional parameters or attributes. Not in use.

- replace:

  Logical; if TRUE, `sds` are recalculated.

- value:

  A vector (or list of vectors) of standard deviations for all stations
  (under each regime and combined).

## Value

[`sds()`](https://tianxia-jia.github.io/mcgf/reference/sds.md) returns
empirical (regime-switching) standard deviations.

## Details

For `mcgf` objects,
[`sds()`](https://tianxia-jia.github.io/mcgf/reference/sds.md) extracts
or computes the empirical standard deviations. The output is a vector of
sds.

For `mcgf_rs` objects,
[`sds()`](https://tianxia-jia.github.io/mcgf/reference/sds.md) extracts
or computes the regime-switching empirical standard deviations. The
output is a list of vectors of sds. Each element in the list corresponds
to the sds for a regime.

`sds<-` assigns `sds` to `x`. Use
[`add_ccfs()`](https://tianxia-jia.github.io/mcgf/reference/ccfs.mcgf.md)
to add both `ccfs` and `sds` to `x`.

## Examples

``` r
data(sim1)
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
sds(sim1_mcgf)
#>         1         2         3         4         5         6         7         8 
#> 0.9862524 0.9748534 0.9770494 0.9698479 0.9599741 0.9994155 0.9646922 0.9475794 
#>         9        10 
#> 0.9896725 0.9755941 

data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
sds(sim2_mcgf)
#> $sds
#>         1         2         3         4         5         6         7         8 
#> 0.9199169 0.9718708 0.9578177 1.0106254 0.9892896 0.9313065 0.9840550 1.0140936 
#>         9        10 
#> 0.9373470 0.9417193 
#> 
#> $sds_rs
#> $sds_rs$`Regime 1`
#>         1         2         3         4         5         6         7         8 
#> 0.9372674 0.9505589 0.9772365 1.0009960 0.9598231 0.9064875 1.0216556 0.9744436 
#>         9        10 
#> 0.9150341 0.9175002 
#> 
#> $sds_rs$`Regime 2`
#>         1         2         3         4         5         6         7         8 
#> 0.8616473 0.9871250 0.9277424 1.0078707 0.9900081 0.9431707 0.9187679 1.0449593 
#>         9        10 
#> 0.9381449 0.9527872 
#> 
#> 
data(sim1)
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim1_sds <- sds(sim1_mcgf)
sds(sim1_mcgf) <- sim1_sds

data(sim2)
sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim2_sds <- sds(sim2_mcgf)
sds(sim2_mcgf) <- sim2_sds
```
