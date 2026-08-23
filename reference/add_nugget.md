# Adjust for nugget effect for correlations

Adjust for nugget effect for correlations

## Usage

``` r
add_nugget(x, nugget)

set_nugget(x, nugget, set_to)
```

## Arguments

- x:

  A correlation matrix or 3-d array.

- nugget:

  A scalar nugget effect.

- set_to:

  A correlation matrix or 3-d array of the same dimension as `x`.

## Value

Correlations of the same dimension as `x`.

## Details

To adjust spatial nugget effect, every entry of `x` is first multiplied
by by \\(1-\text{nugget})\\; Then `add_nugget` adds `nugget` to the
diagonals (or the diagonals of each matrix slice) of `x`, and
`set_nugget` set the diagonals (or the diagonals of each matrix slice)
to the corresponding diagonals of `set_to`.
