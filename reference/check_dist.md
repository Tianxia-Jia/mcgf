# Check if valid distance

Check if valid distance

## Usage

``` r
check_dist(x, name = "x", check_sym = TRUE)
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

Check if `x` is a valid distance vector, matrix or array. It errors if
any elements in `x` is negative, or if `x` is not a symmetric matrix or
an array of symmetric matrices.
