# Calculate (signed) distances between coordinates

Calculate (signed) distances between coordinates

## Usage

``` r
find_dists(locations, longlat = TRUE, origin = 1L, return_grid = FALSE, ...)
```

## Arguments

- locations:

  A matrix or data.frame of 2D points, the first column is x/longitude,
  and the second column is y/latitude.

- longlat:

  Logical, if TRUE Great Circle (WGS84 ellipsoid) distance; if FALSE,
  Euclidean distance.

- origin:

  Optional; used when `longlat` is TRUE. An integer index indicating the
  reference location which will be used as the origin.

- return_grid:

  Logical; used when `longlat` is TRUE. If TRUE the mapped coordinates
  on a 2D plane is returned.

- ...:

  Optional arguments passed to
  [`.find_dists()`](https://tianxia-jia.github.io/mcgf/reference/dot-find_dists.md).

## Value

A list of distance matrices. If `return_grid` is TRUE, a list consists
of a list of distance matrices, the mapped 2D grid, and the origin is
returned.

## Details

`locations` must be a matrix or data.frame containing 2 columns, first
column x/longitude, and second column y/latitude. The row names of
`locations` are used as the names of the locations.

If `longlat` is TRUE, the original coordinates are mapped to a 2D
Euclidean plane given the reference location. First, the Great Circle
(WGS84 ellipsoid) signed distance matrices are calculated, where the
original latitudes are replaced by the the mean of them to find the
signed longitudinal distances and the original longitudes are replaced
by the the mean of them to find the signed latitudinal distances. Then
given the index of a reference location `origin`, a new set of
coordinates in a 2D plane is generated where the coordinates are
determined by the signed distances between the locations and the
reference location. Finally distance matrices of the new coordinates are
outputted.

## Examples

``` r
lon <- c(110, 120, 130)
lat <- c(50, 55, 60)
locations <- cbind(lon, lat)
rownames(locations) <- paste("Site", 1:3)
find_dists(locations)
#> $h
#>           Site 1   Site 2   Site 3
#> Site 1    0.0000 847.5798 1692.986
#> Site 2  847.5798   0.0000  845.408
#> Site 3 1692.9859 845.4080    0.000
#> 
#> $h1
#>           Site 1    Site 2     Site 3
#> Site 1    0.0000 -639.3934 -1275.5045
#> Site 2  639.3934    0.0000  -636.1111
#> Site 3 1275.5045  636.1111     0.0000
#> 
#> $h2
#>          Site 1    Site 2     Site 3
#> Site 1    0.000 -556.3880 -1113.2338
#> Site 2  556.388    0.0000  -556.8459
#> Site 3 1113.234  556.8459     0.0000
#> 
```
