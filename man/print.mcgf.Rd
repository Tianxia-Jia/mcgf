% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print.mcgf.R
\name{print.mcgf}
\alias{print.mcgf}
\title{Print an \code{mcgf} object.}
\usage{
\method{print}{mcgf}(x, attr = ".Data", ...)
}
\arguments{
\item{x}{An \code{mcgf} object.}

\item{attr}{Attribute to be printed.}

\item{...}{Optional arguments to print methods.}
}
\value{
No return value, called for side effects.
}
\description{
Print an \code{mcgf} object.
}
\examples{
data(sim1)
sim1_mcgf <- mcgf(sim1$data, dists = sim1$dists)
print(sim1_mcgf, "dists")
}
