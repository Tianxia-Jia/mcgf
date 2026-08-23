# Simulated regime-switching Markov chain Gaussian field

Simulated RS-MCGF for 20 locations.

## Usage

``` r
sim3
```

## Format

`sim3`: a list containing a data.frame with 5000 rows and 20 columns and
a list of locations.

## Details

`sim3` contains a simulated RS-MCGF for 20 locations. It is simulated
with the same base model and a regime-switching Lagrangian model. The
true parameters for the base model are: \\\text{nugget} = 0, c = 0.05,
\gamma = 0.5, a = 0.5, \alpha = 0.2\\, and the true parameters for the
Lagrangian model are \$\$\text{Regime 1}: \lambda = 0.2, v_1 = -100, v_2
= 100, k = 2,\$\$ \$\$\text{Regime 1}: \lambda = 0.2, v_1 = 200, v_2 =
200, k = 2.\$\$ For parameter estimation, the base model is assumed
known and is used to estimate the regime-switching prevailing winds.

## See also

Other (simulated) datasets:
[`sim1`](https://tianxia-jia.github.io/mcgf/reference/sim1.md),
[`sim2`](https://tianxia-jia.github.io/mcgf/reference/sim2.md),
[`wind`](https://tianxia-jia.github.io/mcgf/reference/wind.md)

## Examples

``` r
# Code used to generate `sim3`
# \donttest{
library(mcgf)
set.seed(123)
x <- stats::rnorm(10, -110)
y <- stats::rnorm(10, 50)
locations <- cbind(x, y)
h <- find_dists(locations, longlat = TRUE)

# simulate regimes
K <- 2
N <- 1000
lag <- 5

tran_mat <- matrix(rnorm(K^2, mean = 0.06 / (K - 1), sd = 0.01), nrow = K)
diag(tran_mat) <- rnorm(K, mean = 0.94, sd = 0.1)
tran_mat <- sweep(abs(tran_mat), 1, rowSums(tran_mat), `/`)
tran_mat
#>            [,1]       [,2]
#> [1,] 0.94635675 0.05364325
#> [2,] 0.06973429 0.93026571
# [,1]       [,2]
# [1,] 0.94635675 0.05364325
# [2,] 0.06973429 0.93026571

regime <- rep(NA, N)
regime[1] <- 1

for (n in 2:(N)) {
    regime[n] <- sample(1:K, 1, prob = tran_mat[regime[n - 1], ])
}
table(regime)
#> regime
#>   1   2 
#> 567 433 
# regime
#   1   2
# 567 433

# simulate RS MCGF
par_base <- list(
    par_s = list(nugget = 0, c = 0.05, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)

par_lagr1 <- list(v1 = -100, v2 = 100, k = 2)
par_lagr2 <- list(v1 = 200, v2 = 200, k = 2)

sim3 <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"), lagrangian_ls = list("lagr_tri"),
    par_base_ls = list(par_base),
    par_lagr_ls = list(par_lagr1, par_lagr2),
    lambda_ls = list(0.2, 0.2),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
sim3 <- sim3[-c(1:(lag + 1)), ]
rownames(sim3) <- 1:nrow(sim3)

sim3 <- list(
    data = sim3[, -1], locations = locations, dists = h,
    label = sim3[, 1]
)
# }
```
