% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcgf.R
\name{new_mcgf}
\alias{new_mcgf}
\title{Create an mcgf object}
\usage{
new_mcgf(data, locations, dists, time, longlat = TRUE, origin = 1L)
}
\arguments{
\item{data}{Time series data set in space-wide format.}

\item{locations}{A matrix of data.frame of 2D points, first column
x/longitude, second column y/latitude. Required when \code{dists} is not supplied.
If longitudes and latitudes are provided, they are mapped to a 2D Euclidean.
See \code{\link[=find_dists]{find_dists()}} for more details.}

\item{dists}{List of signed distance matrices on a 2D Euclidean Plane.
Required when \code{locations} is not supplied.}

\item{time}{Optional, a vector of equally spaced time stamps.}

\item{longlat}{Logical, if TURE \code{locations} contains longitudes and
latitudes.}

\item{origin}{Optional; used when \code{longlat} is TRUE. An integer index
indicating the reference location which well be used as the origin.}
}
\value{
An S3 object of class \code{mcgf}. As it inherits and extends the
\code{data.frame} class, all methods remain valid to the \code{data} part of the
object. Additional attributes may be assigned and extracted.
}
\description{
Create an mcgf object
}
\keyword{internal}
