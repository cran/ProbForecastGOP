\name{matern.variog}
\alias{matern.variog}
\title{Internal function used to estimate the nugget effect, variance, range, and smoothness parameters in the Whittle-Matern variogram. This function is called by the Variog.fit and ProbForecastGOP functions
when fitting a parametric Whittle-Matern variogram to an empirical variogram via the method of Weighted Least Squares.}
\description{
Internal function called by the Variog.fit and the ProbForecastGOP functions when fitting a Whittle-Matern theoretical variogram to an empirical variogram via the method of Weighted Least
Squares.}

\usage{
matern.variog(v,variog,d,w)
}
\arguments{
  \item{v}{numeric vector with the nugget effect, variance, range, and smoothness parameters of the Whittle-Matern variogram.}
  \item{variog}{numeric vector giving the values of the empirical variogram at distances given by the numeric vector d.}
  \item{d}{numeric vector giving the distances (or the bin midpoints) at which the empirical variogram has been computed.}
  \item{w}{numeric vector giving the weights to be used in the Weighted Least Squares, that is, the number of pairs of meteorological stations with distance falling in any given bin.}
}
 
\details{
This function is an internal function that is used and called by the Variog.fit and ProbForecastGOP functions to estimate the parameters of a Whittle-Matern variogram when fitting to an
empirical variogram via Weighted Least Squares.
\bold{- Defaults -}
  
None.
}

\value{
The function returns the weighted least-square loss function relative to the empirical variogram and the theoretical Whittle-Matern variogram evluated for a given set v of nugget effect, variance,
range, and smoothness parameters. This is the function that is minimized to obtain estimates of the nugget effect, variance, range, and smoothness parameters of the Whittle-Matern variogram.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





