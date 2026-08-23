# Check if valid input length

Check if valid input length

## Usage

``` r
check_length_ls(x_ls, length, name)
```

## Arguments

- x_ls:

  List of scaler or vector

- length:

  List of length of `x_ls`.

- name:

  Name of `x` for displaying errors.

## Value

`x_ls`.

## Details

Check if elements in `x_ls` have appropriate length. If length of any
elements in `x_ls` is 1 then they are replicated to match `length`. If
length of any elements is neither 1 or `length` then an error is
signaled.
