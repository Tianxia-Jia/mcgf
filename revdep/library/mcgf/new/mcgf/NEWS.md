## version 1.1.0

---


### Miscellaneous

- New vignettes for regime-switching models
- Updated description and modified reference
- Updated simulated samples throughout the package to allow new locations


### New functions

- R/krige_new.R
- R/find_dists_new.R
- R/cor_lagr_exp.R


### Function updates

- added exponential function for the Lagrangian model
- added two arguments to mcgf()


### Bug fixes

- `find_dists`: fixed bug for computing long/lat distnaces
- `add_lagr`: `dists_base` -> `dists_lagr`


### New features

- Kriging for new locations are supported
- Exponential Lagrangian correlation function added


## version 1.0.1

---


### New functions

- misc/update_NEWS.R


### Function updates

- added arguments to `check_dist_sign` and `check_dist`


### Bug fixes

- `fit_lagr`: fixed check_dist for `dists_lagr`
- `add_lagr`: `dists_base` -> `dists_lagr`


## version 1.0.0

---

- Initial CRAN submimssion

## version 0.0.0.9000

---

### NEWS.md setup

- added NEWS.md creation with [newsmd](https://github.com/Dschaykib/newsmd)

