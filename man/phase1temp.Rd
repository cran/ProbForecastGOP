\name{phase1temp}
\docType{data}
\alias{phase1temp}
\title{Temperature in the Pacific Northwest during January-June 2000}
\description{
 This data set gives the observed and the 48-hour forecast temperature 
(in degree Celsius) at 0 hours GMT for 1092 stations located in the North 
American Pacific Northwest during the period January-June 2000, as 
provided by the Aviation (AVN) member of the University of Washington 
MM5 Mesoscale Ensemble.
}

\usage{data(phase1temp)}

\format{
A matrix with 56489 rows and 6 columns. The columns are:

\emph{date.obs} a numeric vector giving the date of observations in the YYYYMMDDHH format.

\emph{name.stat} a character string giving the identifier of the metereological stations.

\emph{lon.stat} a numeric vector giving the longitudes of the metereological stations.

\emph{lat.stat} a numeric vector giving the latitudes of the metereological stations.

\emph{forecast} a numeric vector giving the 48-hr temperature forecasts as provided by the Aviation (AVN) member of the University of Washington MM5 ensemble.

\emph{obs} a numeric vector giving the observed temperatures.
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
mesoscale weather field forecasting: The Geostatistical Output 
Perturbation (GOP) method (with discussion). \emph{Journal of the American 
Statistical Association}, \bold{Vol. 99 (467),} 575--583.

Grimit, E. P., Mass, C. F. (2002). Initial results of a mesoscale 
short-range ensemble forecasting system over the Pacific Northwest. 
\emph{Weather and Forecasting} \bold{17}, 192--205.}
} 
\keyword{datasets}
