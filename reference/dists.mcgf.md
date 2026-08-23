# Calculating distance matrices for an `mcgf` object

Calculating distance matrices for an `mcgf` object

## Usage

``` r
# S3 method for class 'mcgf'
dists(x, return_grid = FALSE, ...)

dists(x) <- value
```

## Arguments

- x:

  An `mcgf` object.

- return_grid:

  Logical; used when `locations` in `x` are longitudes and latitudes.

- ...:

  Additional parameters or attributes.

- value:

  List of signed distance matrices, outputted from
  [`dists()`](https://tianxia-jia.github.io/mcgf/reference/dists.md).

## Value

A list of signed distance matrices: `h` (Euclidean), `h1` (horizontal),
and `h2` (vertical).

## Details

If the `dists` attribute is available in `x`, it will be printed.
Otherwise `dists` will be calculated based on the `locations` attribute.

## Examples

``` r
data <- cbind(S1 = 1:5, S2 = 4:8, S3 = 5:9)
lon <- c(110, 120, 130)
lat <- c(50, 55, 60)
locations <- cbind(lon, lat)

# if locations are longitudes and latitudes
obj <- mcgf(data = data, locations = locations)
#> `time` is not provided, assuming rows are equally spaced temporally.
obj
#>   S1 S2 S3
#> 1  1  4  5
#> 2  2  5  6
#> 3  3  6  7
#> 4  4  7  8
#> 5  5  8  9
#> 
#> Other attributes: locations, longlat, origin 
dists(obj)
#> $h
#>           1        2        3
#> 1    0.0000 847.5798 1692.986
#> 2  847.5798   0.0000  845.408
#> 3 1692.9859 845.4080    0.000
#> 
#> $h1
#>           1         2          3
#> 1    0.0000 -639.3934 -1275.5045
#> 2  639.3934    0.0000  -636.1111
#> 3 1275.5045  636.1111     0.0000
#> 
#> $h2
#>          1         2          3
#> 1    0.000 -556.3880 -1113.2338
#> 2  556.388    0.0000  -556.8459
#> 3 1113.234  556.8459     0.0000
#> 
dists(obj) <- dists(obj)
obj
#>   S1 S2 S3
#> 1  1  4  5
#> 2  2  5  6
#> 3  3  6  7
#> 4  4  7  8
#> 5  5  8  9
#> 
#> Other attributes: locations, longlat, origin, dists 

# if locations are just coordinates in a 2D plane:
obj <- mcgf(data = data, locations = locations, longlat = FALSE)
#> `time` is not provided, assuming rows are equally spaced temporally.
obj
#>   S1 S2 S3
#> 1  1  4  5
#> 2  2  5  6
#> 3  3  6  7
#> 4  4  7  8
#> 5  5  8  9
#> 
#> Other attributes: locations, longlat, origin 

# calculate distances
dists(obj)
#> $h
#>          1        2        3
#> 1  0.00000 11.18034 22.36068
#> 2 11.18034  0.00000 11.18034
#> 3 22.36068 11.18034  0.00000
#> 
#> $h1
#>    1   2   3
#> 1  0 -10 -20
#> 2 10   0 -10
#> 3 20  10   0
#> 
#> $h2
#>    1  2   3
#> 1  0 -5 -10
#> 2  5  0  -5
#> 3 10  5   0
#> 

# add distances to the `mcgf` object
dists(obj) <- dists(obj)
obj
#>   S1 S2 S3
#> 1  1  4  5
#> 2  2  5  6
#> 3  3  6  7
#> 4  4  7  8
#> 5  5  8  9
#> 
#> Other attributes: locations, longlat, origin, dists 
```
