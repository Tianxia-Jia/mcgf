# Check if valid dists attribute for an `mcgf` object

Check if valid dists attribute for an `mcgf` object

## Usage

``` r
check_dists(dists, n_var, names, name_dists = "dists")
```

## Arguments

- dists:

  List of scaler or vector

- n_var:

  Scalar, number of variables.

- names:

  column and row names of matrices in `dists`.

- name_dists:

  name_dists of `dists`.

## Value

`dists`.

## Details

Check if `dists` is a valid dists attribute for an `mcgf` object. It
errors if 1) `dists` does not contain `h1` or `h2`, 2) if their
dimensions do not match, 3) if it contains elements other than `h`, `h1`
or `h2`. `h` will be computed if it is not given.
