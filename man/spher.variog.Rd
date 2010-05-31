\name{spher.variog}
\alias{spher.variog}
\title{Internal function used to estimate the nugget effect, variance, and range parameters in the spherical variogram. This function is called by the Variog.fit and ProbForecastGOP functions
when fitting a parametric spherical variogram to an empirical variogram via the method of Weighted Least Squares.}
\description{
Internal function called by the Variog.fit and the ProbForecastGOP functions when fitting a spherical theoretical variogram to an empirical variogram via the method of Weighted Least
Squares.}

\usage{
spher.variog(v,variog,d,w)
}
\arguments{
  \item{v}{numeric vector with the nugget effect, variance, and range parameters of the spherical variogram.}
  \item{variog}{numeric vector giving the values of the empirical variogram at distances given by the numeric vector d.}
  \item{d}{numeric vector giving the distances (or the bin midpoints) at which the empirical variogram has been computed.}
  \item{w}{numeric vector giving the weights to be used in the Weighted Least Squares, that is, the number of pairs of meteorological stations with distance falling in any given bin.}
}
 
\details{
This function is an internal function that is used and called by the Variog.fit and ProbForecastGOP functions to estimate the parameters of a spherical variogram when fitting to an
empirical variogram via Weighted Least Squares.
\bold{- Defaults -}
  
None.
}

\value{
The function returns the weighted least-square loss function relative to the empirical variogram and the theoretical spherical variogram evaluated for a given set v of nugget effect, variance, and
range parameters. This is the function that is minimized to obtain estimates of the nugget effect, variance, and range parameters of the spherical variogram.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





