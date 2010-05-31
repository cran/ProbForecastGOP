\name{gencauchy.variog}
\alias{gencauchy.variog}
\title{Internal function used to estimate the nugget effect, variance, range, smoothness, and long-range parameters in the generalized Cauchy variogram. This function is called by the Variog.fit and ProbForecastGOP functions
when fitting a parametric generalized Cauchy variogram to an empirical variogram via the method of Weighted Least Squares.}
\description{
Internal function called by the Variog.fit and ProbForecastGOP functions when fitting a generalized Cauchy theoretical variogram to an empirical variogram via the method of Weighted Least
Squares.}

\usage{
gencauchy.variog(v,variog,d,w)
}
\arguments{
  \item{v}{numeric vector with the nugget effect, variance, range, smoothness, and long-range parameters of the generalized Cauchy variogram.}
  \item{variog}{numeric vector giving the values of the empirical variogram at distances given by the numeric vector d.}
  \item{d}{numeric vector giving the distances (or the bin midpoints) at which the empirical variogram has been computed.}
  \item{w}{numeric vector giving the weights to be used in the Weighted Least Squares, that is, the number of pairs of meteorological stations with distance falling in any given bin.}
}
 
\details{
This function is an internal function that is used and called by the Variog.fit and ProbForecastGOP functions to estimate the parameters of a generalized Cauchy variogram when fitting to an
empirical variogram via Weighted Least Squares.
\bold{- Defaults -}
  
None.
}

\value{
The function returns the weighted least-square loss function relative to the empirical variogram and the theoretical generalized Cauchy variogram evaluated for a given set v of nugget effect, variance,
range, smoothness, and long-range parameters. This is the function that is minimized to obtain estimates of the nugget effect, variance, range, smoothness, and long-range parameters of the generalized Cauchy variogram.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





