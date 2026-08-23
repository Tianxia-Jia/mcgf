# Covariance for joint distribution

Covariance for joint distribution

## Usage

``` r
cov_joint(cov)

cov_par(cov, horizon = 1, n_var, joint = FALSE)
```

## Arguments

- cov:

  Array of covariance matrices.

- horizon:

  Forecast horizon, default is 1.

- n_var:

  Number of locations.

- joint:

  Logical; True if `cov` is the joint covariance matrix.

## Value

The joint covariance matrix for the joint distribution of the current
values and the past values for a Markov chain Gaussian field.

## Details

The covariance matrix of the joint distribution has the block toeplitz
structure. Input `cov` is assumed to be an array of cross-covariance
matrices where the \\i\\th matrix slice correspond to the \\(i-1)\\th
time lag. For example, `cov[, , 1]` is the cross-covariance matrix for
time lag 0. All matrices in `cov` are used to construct the joint
covariance matrix.

`cov_par` gives weights and covariance matrix for the current values.
