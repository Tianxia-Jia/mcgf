# Display fitted models for an `mcgf` or `mcgf_rs` object

Display fitted models for an `mcgf` or `mcgf_rs` object

## Usage

``` r
# S3 method for class 'mcgf'
model(
  x,
  model = c("all", "base", "lagrangian"),
  old = FALSE,
  print_model = TRUE,
  ...
)

# S3 method for class 'mcgf_rs'
model(
  x,
  model = c("all", "base", "lagrangian"),
  old = FALSE,
  print_model = TRUE,
  ...
)
```

## Arguments

- x:

  An mcgf object.

- model:

  Which model to display.

- old:

  Logical; TRUE if the old model needs to be printed.

- print_model:

  Logical; TRUE if time lag and forecast horizon need to be printed.

- ...:

  Additional arguments. Not in use.

## Value

None (invisible `NULL`).

## Details

For `mcgf` and `mcgf_rs` objects,
[`model()`](https://tianxia-jia.github.io/mcgf/reference/model.md)
displays the fitted models and their parameters. When `old = TRUE`, the
old model is printed as well. Note that the old model is not used for
parameter estimation or for kriging.
