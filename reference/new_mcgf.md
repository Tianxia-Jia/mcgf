# Create an mcgf object

Create an mcgf object

## Usage

``` r
new_mcgf(data, locations, dists, time, longlat = TRUE, origin = 1L)
```

## Arguments

- data:

  Time series data set in space-wide format.

- locations:

  A matrix of data.frame of 2D points, first column x/longitude, second
  column y/latitude. Required when `dists` is not supplied. If
  longitudes and latitudes are provided, they are mapped to a 2D
  Euclidean. See
  [`find_dists()`](https://tianxia-jia.github.io/mcgf/reference/find_dists.md)
  for more details.

- dists:

  List of signed distance matrices on a 2D Euclidean Plane. Required
  when `locations` is not supplied.

- time:

  Optional, a vector of equally spaced time stamps.

- longlat:

  Logical, if TRUE `locations` contains longitudes and latitudes.

- origin:

  Optional; used when `longlat` is TRUE. An integer index indicating the
  reference location which will be used as the origin.

## Value

An S3 object of class `mcgf`. As it inherits and extends the
`data.frame` class, all methods remain valid to the `data` part of the
object. Additional attributes may be assigned and extracted.
