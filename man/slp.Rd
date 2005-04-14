\name{slp} 
\docType{data} 
\alias{slp} 
\title{Sea Level Pressure in the Pacific Northwest 
for the period Jan. 12, 2000 - Jan. 14, 2000} 

\description{This data set gives 48-hour forecast and simulated observed 
Sea Level Pressure (measured in millibars, mb) for 59 stations located in 
the Pacific Northwest. 
The 48-hour forecasts of Sea Level Pressure were provided by the Aviation (AVN) 
member of the University of Washington MM5 Mesoscale Ensemble. The observed 
Sea Level Pressure values were simulated for each day using a Multivariate 
Normal distribution with exponential covariance structure. The parameters used 
for the exponential covariance structure were: nugget effect, equal to 
\eqn{0.1} \eqn{mb^2}, \eqn{\sigma^2} equal to \eqn{0.5} 
\eqn{mb^2} and range equal to \eqn{50} km.}
 
\usage{data(slp)} 
\format{A matrix with 150 rows and 6 columns. The columns are:

\emph{date.obs} a numeric vector giving the date of observations in the 
YYYYMMDDHH format.

\emph{id.stat} a character string giving the identifier of the 
metereological stations.

\emph{lon.stat} a numeric vector giving the longitudes of the 
metereological stations.

\emph{lat.stat} a numeric vector giving the latitudes of the 
metereological stations.

\emph{forecast} a numeric vector giving the 48-hr Sea Level Pressure 
forecasts as provided by the Aviation (AVN) member of the University of 
Washington MM5 ensemble.

\emph{obs} a numeric vector giving the simulated Sea Level Pressures.
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated
mesoscale weather field forecasting: The Geostatistical Output
Perturbation (GOP) method (with discussion). 
\emph{Journal of the American Statistical Association}, \bold{Vol. 99 
(467)}, 575--583.

Grimit, E. P., Mass, C. F. (2002). Initial results of a mesoscale 
short-range ensemble forecasting system over the Pacific Northwest. 
\emph{Weather and Forecasting} \bold{17}, 192--205.}
}
\keyword{datasets}

