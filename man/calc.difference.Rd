\name{calc.difference}
\alias{calc.difference}
\title{Internal function used to compute the empirical variogram of a random variable averaged over time.}
\description{
Calculates the square difference between two values of a random variable measured at two meteorological stations. 
This is an internal function used in the estimation of the empirical variogram of a random variable.
}

\usage{
calc.difference(obs)
}
\arguments{
  \item{obs}{numeric vector containing the observations for a random variable at different meteorological stations.}
}
 
\details{
This function is an internal function that is used and called by the avg.variog and the avg.variog.dir functions when
estimating the empirical variogram of a random variable observed at different meteorological stations.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a square matrix of dimension n by n, where n is the length of the vector obs, input to the calc.difference function.
The (i,j) entry of this matrix gives the squared difference between the value of the random variable at meteorological station i, that is, obs_i, and 
the value of the meteorological station j, that is, obs_j.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





