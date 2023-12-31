% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_nugget.R
\name{add_nugget}
\alias{add_nugget}
\alias{set_nugget}
\title{Adjust for nugget effect for correlations}
\usage{
add_nugget(x, nugget)

set_nugget(x, nugget, set_to)
}
\arguments{
\item{x}{A correlation matrix or 3-d array.}

\item{nugget}{A scalar nugget effect.}

\item{set_to}{A correlation matrix or 3-d array of the same dimension as \code{x}.}
}
\value{
Correlations of the same dimension as \code{x}.
}
\description{
Adjust for nugget effect for correlations
}
\details{
To adjust spatial nugget effect, enery entry of \code{x} is first multipled by
by \eqn{(1-\text{nugget})}; Then \code{add_nugget} adds \code{nugget} to the diagonals
(or the diagonals of each matrix slice) of \code{x}, and \code{set_nugget} set the
diagonals (or the diagonals of each matrix slice) to the corresponding
diagonals of \code{set_to}.
}
\keyword{internal}
