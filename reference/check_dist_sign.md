# Check if valid signed distance

Check if valid signed distance

## Usage

``` r
check_dist_sign(x, name, check_sym = TRUE)
```

## Arguments

- x:

  Distance matrix or array.

- name:

  Name of `x` for displaying errors.

- check_sym:

  Logical; if TRUE each matrix (slice) must be symmetric.

## Value

NULL.

## Details

Check if `x` is a valid signed distance vector, matrix or array. It
errors if `x` in absolute value is not a symmetric matrix or an array of
symmetric matrices.
