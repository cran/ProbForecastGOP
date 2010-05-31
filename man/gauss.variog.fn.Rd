\name{gauss.variog.fn}
\alias{gauss.variog.fn}
\title{Internal function used to estimate the variance and range parameters in the gaussian variogram. This function is called by the Variog.fit and ProbForecastGOP functions
when fitting a parametric gaussian variogram to an empirical variogram via the method of Weighted Least Squares 
keeping the nugget effect fixed.}
\description{
Internal function called by the Variog.fit and ProbForecastGOP functions when fitting a gaussian theoretical variogram to an empirical variogram via the method of Weighted Least
Squares when the estimate of the nugget effect is held fixed.}

\usage{
gauss.variog.fn(v,variog,d,w)
}
\arguments{
  \item{v}{numeric vector with the variance and range parameters of the gaussian variogram.}
  \item{variog}{numeric vector giving the values of the empirical variogram at distances given by the numeric vector d.}
  \item{d}{numeric vector giving the distances (or the bin midpoints) at which the empirical variogram has been computed.}
  \item{w}{numeric vector giving the weights to be used in the Weighted Least Squares, that is, the number of pairs of meteorological stations with distance falling in any given bin.}
}
 
\details{
This function is an internal function that is used and called by the Variog.fit and ProbForecastGOP functions to estimate the parameters of a gaussian variogram when fitting to an
empirical variogram via Weighted Least Squares, keeping the nugget effect fixed.
\bold{- Defaults -}
  
None.
}

\value{
The function returns the weighted least-square loss function relative to the empirical variogram and the theoretical gaussian variogram evaluated for a given set v of variance
and range parameters. This is the function that is minimized to obtain estimates of the variance and range parameters of the
gaussian variogram, when the nugget effect is held fixed.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





