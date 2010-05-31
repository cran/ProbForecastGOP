\name{ProbForecastGOP-package}
\Rdversion{1.1}
\alias{ProbForecastGOP-package}
\alias{ProbForecastGOP-package}
\docType{package}
\title{
Probabilistic weather forecast using the GOP method
}
\description{
The ProbForecastGOP package contains a main function, called ProbForecastGOP and other functions, to produce 
probabilistic weather forecasts of weather fields using the Geostatistical Output Perturbation (GOP) method of Gel, Raftery, and Gneiting (JASA, 2004).
}
\details{
\tabular{ll}{
Package: \tab ProbForecastGOP\cr
Type: \tab Package\cr
Version: \tab 1.3.2\cr
Date: \tab 2010-05-31\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

For an overview of how to use the package, including the most important 
functions, please refer to the PDF file describing the package.
}
\author{
Veronica J. Berrocal <veroberrocal@gmail.com>, Yulia Gel, Adrian E. Raftery, Tilmann Gneiting

Maintainer: Veronica J. Berrocal <veroberrocal@gmail.com>
}
\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: The Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Gel, Y., Raftery, A. E., Gneiting, T., Berrocal, V. J. (2004). 
Rejoinder. \emph{Journal of the American Statistical Association}, 
\bold{Vol. 99 (467)}, 588--590.  
}
\keyword{ package }
\seealso{
}
\examples{
library(fields)
library(RandomFields)
data(slp)
day <-slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <-slp$forecast

data(gridlong)
coord1.grid <- gridlong$gridded.lon

data(gridlat)
coord2.grid <- gridlat$gridded.lat

data(forecast.grid)
forecast.grid <- forecast.grid$gridded.forecast

## Specified cutpoints, default values for all the other fields.
## Only empirical variogram computation
empirical.variog <- ProbForecastGOP(day=day,obs=obs,forecast=forecast,
id=id,coord1=coord1,coord2=coord2,cut.points=seq(0,1000,by=2),
max.dist=NULL,nbins=NULL,variog.model="exponential",max.dist.fit=NULL,
init.val=NULL,fix.nugget=FALSE,coord1.grid=coord1.grid,
coord2.grid=coord2.grid,forecast.grid=forecast.grid,n.sim=10,
out="VARIOG",n.displ=4,qt.displ=c(5,50,95))


## Unspecified cutpoints.
## Fit of a parametric variogram to an empirical variogram.
fitted.variog <- ProbForecastGOP(day=day,obs=obs,forecast=forecast,id=id,
coord1=coord1,coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL,
variog.model="exponential",max.dist.fit=NULL,init.val=NULL,fix.nugget=FALSE,
coord1.grid=coord1.grid,coord2.grid=coord2.grid,forecast.grid=forecast.grid,
n.sim=10,out="FIT",n.displ=4,qt.displ=c(5,50,95))


## Unspecified cutpoints. 
## Whole routine.
simulation.fields <- 
ProbForecastGOP(day=day,obs=obs,forecast=forecast,id=id,coord1=coord1,
coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL,variog.model=NULL,
max.dist.fit=NULL,init.val=NULL,fix.nugget=FALSE,coord1.grid=coord1.grid,
coord2.grid=coord2.grid,forecast.grid=forecast.grid,n.sim=4,out="SIM",
n.displ=4,qt.displ=c(5,50,95))
}
