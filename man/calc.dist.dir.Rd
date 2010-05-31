\name{calc.dist.dir}
\alias{calc.dist.dir}
\title{Internal function used to compute the inter-site directional distance between sites.}
\description{
Calculates the directional distance between sites of a random variable. This is an internal function called by the avg.variog.dir and the EmpDir.variog function
when computing the directional empirical variogram.
}

\usage{
calc.dist.dir(coord1,coord2,id,tol.angle.rad1,tol.angle.rad2,type)
}
\arguments{
  \item{coord1}{numeric vector containing the longitudes of the meteorological stations.}
  \item{coord2}{vector containing the latitudes of the metereological stations.}
  \item{id}{vector with the id of the metereological stations.}
  \item{tol.angle1}{number giving a lower bound for the tolerance angle (measured in radians).}
  \item{tol.angle2}{number giving an upper bound for the tolerance angle (measured in radians).}
  \item{type}{character string indicating the direction to use for variogram computations. Possible values are either 'E' (for East-West) or 'N' (for North-South).} 
}
 
\details{
This function is an internal function that is used and called by the avg.variog.dir functions when
estimating the directional empirical variogram of a random variable observed at different meteorological stations.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a square matrix of dimension n by n, where n is the length of the vector coord1 where the
(i,j) entry of represents the directional distance between the i-th and the j-th metereological stations.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





