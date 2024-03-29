% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_dists_new.R
\name{.find_dists_new}
\alias{.find_dists_new}
\title{Calculate (signed) distances between coordinates}
\usage{
.find_dists_new(grid, grid_new, names = NULL, names_new = NULL, longlat = TRUE)
}
\arguments{
\item{grid}{A matrix of 2D points, first column x/longitude, second column
y/latitude.}

\item{grid_new}{A matrix of 2D points, first column x/longitude, second column
y/latitude.}

\item{names_new}{Names of new locations.}

\item{longlat}{Logical, if TURE Great Circle (WGS84 ellipsoid) distance;
if FALSE, Euclidean distance.}
}
\value{
List of signed distances between the new locations and the old grid.
}
\description{
Calculate (signed) distances between coordinates
}
\keyword{internal}
