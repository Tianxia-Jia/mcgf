# Create an mcgf_rs object

Create an mcgf_rs object

## Usage

``` r
new_mcgf_rs(x, label)
```

## Arguments

- x:

  An mcgf object.

- label:

  A vector of regime labels. Its length must be the same as the number
  of rows in `data`.

## Value

An S3 object of class `mcgf_rs`. As it inherits and extends the `mcgf`
and then the `data.frame` class, all methods remain valid to the `data`
part of the object. Additional attributes may be assigned and extracted.
