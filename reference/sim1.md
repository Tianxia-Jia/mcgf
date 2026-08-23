# Simulated Markov chain Gaussian field

Simulated MCGF for 10 locations.

## Usage

``` r
sim1
```

## Format

`sim1`: a list containing a data.frame with 1000 rows and 10 columns and
a list of distances

## Details

`sim1` contains a simulated MCGF for 10 locations. It is simulated with
a separable base model and a triangular Lagrangian model. The true
parameters for the base model are: \\\text{nugget} = 0, c = 0.001,
\gamma = 0.5, a = 0.5, \alpha = 0.8\\, and those for the Lagrangian
model are: \\v1 = 200, v2 = 200, k = 2, \lambda = 0.2\\

## See also

Other (simulated) datasets:
[`sim2`](https://tianxia-jia.github.io/mcgf/reference/sim2.md),
[`sim3`](https://tianxia-jia.github.io/mcgf/reference/sim3.md),
[`wind`](https://tianxia-jia.github.io/mcgf/reference/wind.md)

## Examples

``` r
# Code used to generate `sim1`
# \donttest{
library(mcgf)
set.seed(123)
x <- stats::rnorm(10, -110)
y <- stats::rnorm(10, 50)
locations <- cbind(x, y)
h <- find_dists(locations, longlat = TRUE)

N <- 1000
lag <- 5

par_base <- list(
    par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.8)
)
par_lagr <- list(v1 = 200, v2 = 200, k = 2)

sim1 <- mcgf_sim(
    N = N, base = "sep", lagrangian = "lagr_tri",
    par_base = par_base, par_lagr = par_lagr, lambda = 0.2,
    dists = h, lag = lag
)
sim1 <- sim1[-c(1:(lag + 1)), ]
rownames(sim1) <- 1:nrow(sim1)

sim1 <- list(data = sim1, locations = locations, dists = h)
# }
```
