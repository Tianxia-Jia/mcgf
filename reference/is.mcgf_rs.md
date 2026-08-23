# Check if an object is an `mcgf_rs` object.

Check if an object is an `mcgf_rs` object.

## Usage

``` r
is.mcgf_rs(x)

as.mcgf_rs(x, label, ncores = 1)
```

## Arguments

- x:

  An Object.

- label:

  A vector of regime labels. Its length must be the same as the number
  of rows in `data`.

- ncores:

  Number of cpu cores used for computing in `[ccfs()]`.

## Value

`is.mcgf_rs` returns a logical value; TRUE if `x` is of the `mcgf_rs`
class. `as.mcgf_rs` coerces an `mcgf` object to an `mcgf_rs` object by
adding regime labels. Fitted base or Lagrangian models in `x` are kept.

## Examples

``` r
data(sim2)
is.mcgf_rs(sim2)
#> [1] FALSE

sim2_mcgf <- mcgf(sim2$data, dists = sim2$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
is.mcgf_rs(sim2_mcgf)
#> [1] FALSE

sim2_mcgf <- mcgf_rs(sim2$data, dists = sim2$dists, label = sim2$label)
#> `time` is not provided, assuming rows are equally spaced temporally.
is.mcgf_rs(sim2_mcgf)
#> [1] TRUE
data(sim2)
sim2_mcgf <- mcgf(sim2$data, dists = sim2$dists)
#> `time` is not provided, assuming rows are equally spaced temporally.
sim2_mcgf <- as.mcgf_rs(sim2_mcgf, label = sim2$label)
```
