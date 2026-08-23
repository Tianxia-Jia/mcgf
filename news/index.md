# Changelog

## mcgf (development version)

## mcgf 1.2.0

CRAN release: 2026-08-21

### version 1.2.0

------------------------------------------------------------------------

#### Miscellaneous

- Updated package documentation and vignettes for consistency and
  clarity
- Added the reference to Jia and Sezer (2025) for regime-switching
  spatio-temporal covariance models
- Corrected typographical, grammatical, and terminology errors
  throughout the package

#### Function updates

- Standardized support and documentation for the exponential Lagrangian
  correlation function
- Updated documentation for weighted least squares and approximate
  conditional maximum likelihood estimation
- Improved documentation and argument descriptions for covariance,
  fitting, and Kriging functions

#### Bug fixes

- `mcgf_sim`: fixed handling of the default Lagrangian model argument
- Improved covariance-matrix validation for maximum likelihood
  estimation
- Fixed issues in distance assignment and Kriging functions for existing
  and new locations
- Improved input validation and error handling across model fitting and
  forecasting functions

### version 1.1.1

------------------------------------------------------------------------

#### Function updates

- added `dists_new_base` argument to `krige_new.mcgf`

#### Bug fixes

- `.find_dists`: fixed an argument bug for `names`
- `obj_mle`: NA values are removed in llike

### version 1.1.0

------------------------------------------------------------------------

#### Miscellaneous

- New vignettes for regime-switching models
- Updated description and modified reference
- Updated simulated samples throughout the package to allow new
  locations

#### New functions

- R/krige_new.R
- R/find_dists_new.R
- R/cor_lagr_exp.R

#### Function updates

- added exponential function for the Lagrangian model
- added two arguments to mcgf()

#### Bug fixes

- `find_dists`: fixed bug for computing long/lat distnaces
- `add_lagr`: `dists_base` -\> `dists_lagr`

#### New features

- Kriging for new locations are supported
- Exponential Lagrangian correlation function added

### version 1.0.1

------------------------------------------------------------------------

#### New functions

- misc/update_NEWS.R

#### Function updates

- added arguments to `check_dist_sign` and `check_dist`

#### Bug fixes

- `fit_lagr`: fixed check_dist for `dists_lagr`
- `add_lagr`: `dists_base` -\> `dists_lagr`

### version 1.0.0

------------------------------------------------------------------------

- Initial CRAN submimssion

### version 0.0.0.9000

------------------------------------------------------------------------

#### NEWS.md setup

- added NEWS.md creation with
  [newsmd](https://github.com/Dschaykib/newsmd)
