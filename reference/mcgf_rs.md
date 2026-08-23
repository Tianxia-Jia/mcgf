# Create mcgf_rs object

Create mcgf_rs object

## Usage

``` r
mcgf_rs(data, locations, dists, label, time, longlat = TRUE, origin = 1L)
```

## Arguments

- data:

  Time series data set in space-wide format.

- locations:

  A matrix of data.frame of 2D points, first column longitude, second
  column latitude, both in decimal degrees. Required when `dists` is not
  supplied.

- dists:

  List of signed distance matrices. Required when `locations` is not
  supplied.

- label:

  A vector of regime labels. Its length must be the same as the number
  of rows in `data`.

- time:

  Optional, a vector of equally spaced time stamps.

- longlat:

  Logical, if `locations` contains longitudes and latitudes.

- origin:

  Optional; used when `longlat` is TRUE. An integer index indicating the
  reference location which will be used as the origin.

## Value

An S3 object of class `mcgf_rs`. As it inherits and extends the `mcgf`
and then the `data.frame` class, all methods remain valid to the `data`
part of the object. Additional attributes may be assigned and extracted.

## Details

An `mcgf_rs` object extends the S3 classes `mcgf` and `data.frame`.

For inputs, `data` must be in space-wide format where rows correspond to
different time stamps and columns refer to spatial locations. Supply
either `locations` or `dists`. `locations` is a matrix or data.frame of
2D points with first column x/longitude and second column y/latitude. By
default it is treated as a matrix of Earth's coordinates in decimal
degrees. Number of rows in `locations` must be the same as the number of
columns of `data`. `dists` must be a list of signed distance matrices
with names `h1`, `h2`, and `h`. If `h` is not given, it will be
calculated as the Euclidean distance of `h1` and `h2`. `time` is a
vector of equally spaced time stamps. If it is not supplied then `data`
is assumed to be temporally equally spaced. `label` must be a vector
containing regime labels, and its length must be the same as the number
of rows in `x`.

An `mcgf_rs` object extends the S3 classes `mcgf` and `data.frame`, all
methods remain valid to the `data` part of the object.

## Examples

``` r
data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
lon <- c(110, 120, 130)
lat <- c(50, 55, 60)
locations <- cbind(lon, lat)
label <- c(1, 1, 2, 2, 2)
obj <- mcgf_rs(data, locations = locations, label = label)
#> `time` is not provided, assuming rows are equally spaced temporally.
print(obj, "locations")
#>      lon lat
#> [1,] 110  50
#> [2,] 120  55
#> [3,] 130  60
print(obj, "label")
#> [1] 1 1 2 2 2
#> Levels: 1 2
```
