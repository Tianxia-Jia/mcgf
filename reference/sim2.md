# Simulated regime-switching Markov chain Gaussian field

Simulated RS-MCGF for 10 locations.

## Usage

``` r
sim2
```

## Format

`sim2`: a list containing a data.frame with 1000 rows and 10 columns, a
list of distances, and a vector of regime labels.

## Details

`sim2` contains a simulated RS-MCGF for 10 locations. It is simulated
with a regime-switching separable base model. The true parameters for
the base model are: \$\$\text{Regime 1}: \text{nugget} = 0, c = 0.01,
\gamma = 0.5, a = 0.5, \alpha = 0.2,\$\$ \$\$\text{Regime 2}:
\text{nugget} = 0, c = 0.04, \gamma = 0.5, a = 0.3, \alpha = 0.9.\$\$

## See also

Other (simulated) datasets:
[`sim1`](https://tianxia-jia.github.io/mcgf/reference/sim1.md),
[`sim3`](https://tianxia-jia.github.io/mcgf/reference/sim3.md),
[`wind`](https://tianxia-jia.github.io/mcgf/reference/wind.md)

## Examples

``` r
# Code used to generate `sim2`
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
#            [,1]       [,2]
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
par_base1 <- list(
    par_s = list(nugget = 0, c = 0.001, gamma = 0.5),
    par_t = list(a = 0.5, alpha = 0.2)
)

par_base2 <- list(
    par_s = list(nugget = 0, c = 0.004, gamma = 0.5),
    par_t = list(a = 0.3, alpha = 0.9)
)

sim2 <- mcgf_rs_sim(
    N = N, label = regime,
    base_ls = list("sep"), lagrangian_ls = list("none"),
    par_base_ls = list(par_base1, par_base2),
    lambda_ls = list(0.1, 0.3),
    lag_ls = list(lag, lag),
    dists_ls = list(h, h)
)
sim2 <- sim2[-c(1:(lag + 1)), ]
rownames(sim2) <- 1:nrow(sim2)

sim2 <- list(
    data = sim2[, -1], locations = locations, dists = h,
    label = sim2[, 1]
)
# }
```
