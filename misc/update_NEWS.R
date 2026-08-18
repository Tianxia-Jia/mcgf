library(newsmd)

# update roxygen
roxygen2::roxygenise()

# Remove default NEWS.md
unlink("NEWS.md")

#------------------------------------------------------------------------------#
# Initial submission
#------------------------------------------------------------------------------#

my_news <- news$new()
my_news <- newsmd()

#------------------------------------------------------------------------------#
# version 1.0.0
#------------------------------------------------------------------------------#

my_news$add_version("1.0.0")

my_news$add_bullet("Initial CRAN submimssion")

my_news$get_text()
my_news$write()

#------------------------------------------------------------------------------#
# version 1.0.1
#------------------------------------------------------------------------------#

my_news$add_version("1.0.1")

my_news$add_subtitle("Bug fixes")
my_news$add_bullet(c("`fit_lagr`: fixed check_dist for `dists_lagr`",
                     "`add_lagr`: `dists_base` -> `dists_lagr`"))

my_news$add_subtitle("Function updates")
my_news$add_bullet(c("added arguments to `check_dist_sign` and `check_dist`"))

my_news$add_subtitle("New functions")
my_news$add_bullet(c("misc/update_NEWS.R"))

my_news$get_text()
my_news$write()

#------------------------------------------------------------------------------#
# version 1.1.0
#------------------------------------------------------------------------------#

my_news$add_version("1.1.0")

my_news$add_subtitle("New features")
my_news$add_bullet(c("Kriging for new locations are supported",
                     "Exponential Lagrangian correlation function added"))

my_news$add_subtitle("Bug fixes")
my_news$add_bullet(c("`find_dists`: fixed bug for computing long/lat distnaces",
                     "`add_lagr`: `dists_base` -> `dists_lagr`"))

my_news$add_subtitle("Function updates")
my_news$add_bullet(c("added two arguments to mcgf()"))
my_news$add_bullet(c("added exponential function for the Lagrangian model"))

my_news$add_subtitle("New functions")
my_news$add_bullet(c("R/cor_lagr_exp.R"))
my_news$add_bullet(c("R/find_dists_new.R"))
my_news$add_bullet(c("R/krige_new.R"))

my_news$add_subtitle("Miscellaneous")
my_news$add_bullet(c(
    "Updated simulated samples throughout the package to allow new locations"
))
my_news$add_bullet(c("Updated description and modified reference"))
my_news$add_bullet(c("New vignettes for regime-switching models"))

my_news$get_text()
my_news$write()

#------------------------------------------------------------------------------#
# version 1.1.1
#------------------------------------------------------------------------------#

my_news$add_version("1.1.1")

my_news$add_subtitle("Bug fixes")
my_news$add_bullet(c("`obj_mle`: NA values are removed in llike"))
my_news$add_bullet(c("`.find_dists`: fixed an argument bug for `names`"))

my_news$add_subtitle("Function updates")
my_news$add_bullet(c("added `dists_new_base` argument to `krige_new.mcgf`"))

my_news$get_text()
my_news$write()

#------------------------------------------------------------------------------#
# version 1.2.0
#------------------------------------------------------------------------------#

my_news$add_version("1.2.0")

my_news$add_subtitle("Bug fixes")
my_news$add_bullet(c(
    "`mcgf_sim`: fixed handling of the default Lagrangian model argument",
    "Improved covariance-matrix validation for maximum likelihood estimation",
    "Fixed issues in distance assignment and Kriging functions for existing and new locations",
    "Improved input validation and error handling across model fitting and forecasting functions"
))

my_news$add_subtitle("Function updates")
my_news$add_bullet(c(
    "Standardized support and documentation for the exponential Lagrangian correlation function",
    "Updated documentation for weighted least squares and approximate conditional maximum likelihood estimation",
    "Improved documentation and argument descriptions for covariance, fitting, and Kriging functions"
))

my_news$add_subtitle("Testing and documentation")
my_news$add_bullet(c(
    "Expanded `testthat` coverage across object construction, distance calculations, correlation functions, simulation, model fitting, covariance construction, and Kriging",
    "Added user-oriented vignettes for covariance-model selection, simulation, and forecasting at new locations",
    "Updated the README with a hands-on workflow for fitting MCGF and RS-MCGF models",
    "Improved existing vignettes and examples for consistency with the current package interface"
))

my_news$add_subtitle("Miscellaneous")
my_news$add_bullet(c(
    "Added the reference to Jia and Sezer (2025) for regime-switching spatio-temporal covariance models",
    "Corrected typographical, grammatical, and terminology errors throughout the package"
))

my_news$get_text()
my_news$write()

my_news$get_text()
my_news$write()
