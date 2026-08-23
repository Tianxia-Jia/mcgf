# Check if an object is an `mcgf` object.

Check if an object is an `mcgf` object.

## Usage

``` r
is.mcgf(x)
```

## Arguments

- x:

  An Object.

## Value

Logical; TRUE if `x` is of the `mcgf` class

## Examples

``` r
data(sim1)
is.mcgf(sim1)
#> [1] FALSE

sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
is.mcgf(sim1_mcgf)
#> [1] TRUE
```
