# Check if valid input length

Check if valid input length

## Usage

``` r
check_length(x, length, name)
```

## Arguments

- x:

  Scalar or vector

- length:

  Length of `x`.

- name:

  Name of `x` for displaying errors.

## Value

`x`.

## Details

Check if `x` has appropriate length. If length of `x` is 1 then `x` is
replicated to match `length`. If length of `x` is neither 1 nor `length`
then an error is signaled.
