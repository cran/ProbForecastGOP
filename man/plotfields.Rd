\name{plotfields}
\alias{plotfields}
\title{Plot of weather fields}
\description{
Plot weather fields.}

\usage{
plotfields(field, x.lim, y.lim, title)
}
\arguments{
  \item{field}{numeric square matrix containing the values of the weather field on a grid.}
  \item{x.lim}{numeric vector giving the smallest and the largest longitude to be displayed.}
  \item{y.lim}{numeric vector giving the smallest and the largest latitude to be displayed.}
  \item{title}{character string with the title for the plot.}
}

\value{
The function returns a graphical display of the weather field on a region delimited by the lower and upper bound for the longitude and the latitude.
}
\author{Gel, Y., Raftery, A. E., Gneiting, T., Berrocal, V. J. \email{veronica@stat.washington.edu}.}

\note{This function uses the package \pkg{fields}.}

\seealso{The package \pkg{fields} for display of spatial data, and \code{\link[fields]{US}} and \code{\link[fields]{world}} for a map of the US and the world.}

\reference{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: The Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Nychka, D. (2004). \emph{The fields package}. Available at: http:lib.stat.cmu.edu/R/CRAN/doc/package/fields.pdf.}

\examples{
## Loading data
library(fields)
library(RandomFields)
data(slp)
data(gridlong)
data(gridlat)
data(forecast.grid)
day <- slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <- slp$forecast
coord1.grid <- gridlong$gridded.long
coord2.grid <- gridlat$gridded.lat
forecast.grid <- forecast.grid$gridded.forecast

## Computing the empirical variogram
variogram <- Emp.variog(day,obs,forecast,id,coord1,coord2,cut.points=NULL,
max.dist=NULL,nbins=NULL)

## Estimating parameters
param.est <- Variog.fit(variogram,"exponential",max.dist.fit=NULL,
init.val=NULL,fix.nugget=FALSE)

## Simulating realizations of the weather random field
simul <- Field.sim(obs, forecast, coord1.grid, coord2.grid, forecast.grid,
variog.model="exponential", param.est=c(param.est$nugget,param.est$variance,
param.est$range), n.sim=4, n.displ=0, qt.displ=c(10,50,90))

##Plotting one of the simulated weather random fields
par(mfrow=c(1,1))
plotfields(simul$sim.fields[,,1],x.lim=c(min(coord1.grid),max(coord1.grid)),
y.lim=c(min(coord2.grid),max(coord2.grid)),title="Simulated weather field")

## Plotting one of the percentiles of the weather field
par(mfrow=c(1,1))
plotfields(simul$pct.fields[,,1],x.lim=c(min(coord1.grid),max(coord1.grid)),
y.lim=c(min(coord2.grid),max(coord2.grid)),title="10th percentile")

\dontshow{
## The function is currently defined as
function(field,x.lim,y.lim,title){

# Input check
dim.f <- dim(field)

if(length(dim.f)==0){
  stop("Invalid input for field")
  }

if(length(dim.f)!=2){
  stop("Invalid input for field")
  }

if(dim.f[1]!=65 | dim.f[2]!=65){
  stop("Invalid input for field")
  }

dim.f2 <- dim.f^2

if(sum(is.numeric(field))==0){
  stop("Invalid input for field")
  }

## input check for the x vector
l.limx <- length(x.lim)
if(l.limx < 2){
  stop("Invalid input for x.lim")
  }

if(l.limx > 2){
  stop("Invalid input for x.lim")
  }

if(l.limx==2){
  if(sum(is.numeric(x.lim)==rep(TRUE,2))<2){
  stop("Invalid input for x.lim")
  }
  if(x.lim[1] >= x.lim[2]){ 
  stop("Invalid input for x.lim")
  }
  }


## input check on the lims for the y vector 
l.limy <- length(y.lim)
if(l.limy < 2){
  stop("Invalid input for y.lim")
  }

if(l.limy > 2){
  stop("Invalid input for y.lim")
  }

if(l.limy==2){
  if(sum(is.numeric(y.lim)==rep(TRUE,2))<2){
  stop("Invalid input for y.lim")
  }
  if(y.lim[1] >= y.lim[2]){ 
  stop("Invalid input for y.lim")
  }
  }

lims <- c(min(field,na.rm=TRUE),max(field,na.rm=TRUE))
size=65  
hor.crd <- seq(x.lim[1],x.lim[2],length = size)  
ver.crd <- seq(y.lim[1],y.lim[2],length = size)  
n.level <- 5
US.map <- 0
par(ask=FALSE)
image.plot(hor.crd, ver.crd, field, xlim=c(min(hor.crd), max(hor.crd)), ylim = c(min(ver.crd),max(ver.crd)),zlim=lims,legend.width=.015,legend.shrink=.8,
main=title,xlab="",ylab="",offset = 0.02, col=rainbow(100,start=0,end=0.85)[100:1])

if(min(hor.crd) <= -124.7 & max(hor.crd) >= -124.7){
  US.map <- 1
  }

if(min(hor.crd) >= -124.7 & max(hor.crd) <= -67.1){
  US.map <- 1
  }

if(min(hor.crd) <= -67.1 & max(hor.crd) >= -67.1){
  US.map <- 1
  }

if(min(ver.crd) <= 25.2 & max(ver.crd) >= 25.2){
  US.map <- 1 + US.map
  }

if(min(ver.crd) >= 25.2 & max(ver.crd) <= 49.4){
  US.map <- 1 + US.map
  }

if(min(ver.crd) <= 49.4 & max(ver.crd) >= 49.4){
  US.map <- 1 + US.map
  }


if(US.map==2){
  US(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
  world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)     
  }  
   
if(US.map < 2){
  world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
  }

  }
}
}
\keyword{file}
