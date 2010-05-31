\name{calc.dist}
\alias{calc.dist}
\title{Internal function used to compute the inter-site distance between sites.}
\description{
Calculates the distance between sites of a random variable. This is an internal function called by the avg.variog and the Emp.variog function
when computing the empirical variogram.
}

\usage{
calc.dist(coord1,coord2,id)
}
\arguments{
  \item{coord1}{numeric vector containing the longitudes of the meteorological stations.}
  \item{coord2}{vector containing the latitudes of the metereological stations.}
  \item{id}{vector with the id of the metereological stations.}
}
 
\details{
This function is an internal function that is used and called by the avg.variog.dir functions when
estimating the empirical variogram of a random variable observed at different meteorlogocal stations.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a square matrix of dimension n by n, where n is the length of the vector coord1 where the
(i,j) entry of represents the distance between the i-th and the j-th metereological stations.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





