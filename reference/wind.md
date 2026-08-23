# Ireland wind data, 1961-1978

Daily average wind speeds for 1961-1978 at 11 synoptic meteorological
stations in the Republic of Ireland (Haslett and raftery 1989). Wind
speeds are in m/s. De-trended data sets are also provided.

## Usage

``` r
wind
```

## Format

`wind`: a list containing a data.frame with 6574 rows and 12 columns,
and a list of locations.

## Details

The data were obtained from the **gstat** package, and were modified so
that the first column is the time stamps. Locations of the 11 stations
are given in `wind_loc`. `wind_train` and `wind_test` contain de-trended
and square-root transformed train (1961-1970) and test (1971-1978) data
sets. See Gneiting et al. (2006) for de-trending details. `wind_trend`
contains the estimated annual trend and station-wise mean from the
training dataset.

## References

Haslett, J. and Raftery, A. E. (1989). Space-time Modelling with
Long-memory Dependence: Assessing Ireland's Wind Power Resource (with
Discussion). Applied Statistics 38, 1-50.

Gneiting, T., Genton, M., & Guttorp, P. (2006). Geostatistical
Space-Time Models, Stationarity, Separability, and Full Symmetry. In
C&H/CRC Monographs on Statistics & Applied Probability (pp. 151–175).
Chapman and Hall/CRC.

## See also

Other (simulated) datasets:
[`sim1`](https://tianxia-jia.github.io/mcgf/reference/sim1.md),
[`sim2`](https://tianxia-jia.github.io/mcgf/reference/sim2.md),
[`sim3`](https://tianxia-jia.github.io/mcgf/reference/sim3.md)
