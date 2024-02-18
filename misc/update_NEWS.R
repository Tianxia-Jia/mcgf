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
