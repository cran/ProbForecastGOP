\name{sim.field}
\alias{sim.field}
\title{Internal function called by the Field.sim and the ProbForecastGOP functions to simulate realizations of forecast weather fields.}
\description{
Internal function called by the Field.sim and the ProbForecastGOP functions to simulate realizations of forecast weather fields. It returns a matrix with values of the realizations.}

\usage{
sim.field(variog.model,param.est,x,y,n.sim)
}
\arguments{
  \item{variog.model}{character string with the name of the variogram model to be used for the simulations. Implemented models are \kbd{exponential}, \kbd{spherical}, \kbd{gauss}, \kbd{matern}, and \kbd{gencauchy}.}
  \item{param.est}{numeric vector containing values for the parameters to use in the parametric variogram model.}
  \item{x}{numeric vector of length n containing the longitudes of the grid points at which the the realizations of the weather forecast fields should be given.}
  \item{y}{numeric vector of length n containing the latitudes of the grid points at which the the realizations of the weather forecast fields should be given.}
  \item{n.sim}{number of realizations to be simulated.}
}
 
\details{
This function is an internal function that is used and called by the Field.sim and the ProbForecastGOP function to simulate realizations of forecast weather fields.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a numeric matrix of dimension n by n.sim where each column represents a realization of a forecast weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





