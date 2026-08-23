# Calculate (signed) distances between coordinates

Calculate (signed) distances between coordinates

## Usage

``` r
.find_dists_new(
  grid,
  grid_new,
  names = NULL,
  names_new = NULL,
  longlat = TRUE,
  origin = 1L,
  return_grid = FALSE,
  lon_ref,
  lat_ref
)
```

## Arguments

- grid:

  A matrix of 2D points, first column x/longitude, second column
  y/latitude.

- grid_new:

  A matrix of 2D points, first column x/longitude, second column
  y/latitude.

- names:

  Names of locations.

- names_new:

  Names of new locations.

- longlat:

  Logical, if TRUE Great Circle (WGS84 ellipsoid) distance; if FALSE,
  Euclidean distance.

- origin:

  Optional; used when `longlat` is TRUE. An integer index indicating the
  reference location which will be used as the origin.

- return_grid:

  Logical; used when `longlat` is TRUE. If TRUE the mapped coordinates
  on a 2D plane is returned.

- lon_ref:

  Reference longitude when computing the longitudinal distances. Default
  is the mean of longitudes in `grid`.

- lat_ref:

  Reference latitude when computing the latitudinal distances. Default
  is the mean of latitudes in `grid`.

## Value

List of signed distances between the new locations and the old grid.
