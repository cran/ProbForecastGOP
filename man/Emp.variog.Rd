\name{Emp.variog}
\alias{Emp.variog}
\title{Empirical variogram averaged over time}
\description{
Calculates the empirical variogram of forecast errors, averaged over time. 
}

\usage{
Emp.variog(day, obs, forecast, id, coord1, coord2, 
cut.points=NULL, max.dist=NULL, nbins=300)
}

\arguments{
  \item{day}{numeric vector containing the day of observation.}
  \item{obs}{numeric vector containing the observed weather quantity.}
  \item{forecast}{numeric vector containing the forecasted weather quantity.}
  \item{id}{vector with the id of the metereological stations in the dataset.}
  \item{coord1}{vector containing the longitudes of the metereological stations.}
  \item{coord2}{vector containing the latitudes of the metereological stations.}
  \item{cut.points}{numeric vector containing the cutpoints used for variogram binning.}
  \item{max.dist}{a numerical value giving the upper bound on the distance considered in the variogram computation.}
  \item{nbins}{a numerical value giving the number of bins for variogram binning. If both cut.points and nbins are entered, the entry for nbins will be ignored and the vector with the cutpoints will instead be used for variogram binning.} 
}

\details{
The function includes bias-correction; it regresses the forecast on the 
observed weather quantity and computes the residuals. The empirical variogram of the residuals is then calculated by determining, for each day, the distance among all pairs of stations that have been observed in the same day and by calculating for each day the sum of all the squared differences in the residuals within each bin. These sums are then averaged over time, with weights for each bin given by the sum over time of the number of pairs of stations within the bin. 

The formula used is:
 \deqn{ \gamma(h) = \sum_d \frac{1}{2N_{(h,d)}} (\sum_i (Y(x_{i}+h,d)-Y(x_{i},d))^2)}

where \eqn{\gamma(h)} is the empirical variogram at distance \eqn{h}, \eqn{N_{(h,d)}} is the number of pairs of stations that have been recorded at day \eqn{d} and whose 
distance is equal to \eqn{h}, and \eqn{Y(x_{i}+h,d)} and \eqn{Y(x_{i},d)} 
are, respectively, the value of weather quantity observed on day \eqn{d} at stations located at \eqn{x_{i}+h} and \eqn{x_{i}}. Variogram binning is ignored in this formula. 

- Defaults -
 
If the vector with the cutpoints is not specified, the cutpoints are determined so that there are \code{nbins} bins with approximately the same number of pairs per bin. 

If both the vector with the cutpoints and the number of bins, \code{nbins}, are unspecified, the function by default determines the cutpoints so that there are 300 bins with approximately the same number of pairs per bin. If both the vector with the cutpoints and the number of bins are provided, the entry for the number of bins is ignored and the vector with the cutpoints is used for variogram binning.

The default value for the maximum distance considered in the variogram 
computation is the \eqn{90}-th percentile of the distances between the stations. 
}

\value{
The function returns a list with components given by:
  \item{res.var}{Variance of the forecast errors.}
  \item{bin.midpoints}{Numeric vector with midpoints of the bins used in the empirical variogram computation.}
  \item{number.pairs}{Numeric vector with the number of pairs per bin.}
  \item{empir.variog}{Numeric vector with the empirical variogram values.}
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated probabilistic 
mesoscale weather field forecasting: The Geostatistical Output 
Perturbation (GOP) method (with discussion). \emph{Journal of the American 
Statistical Association}, \bold{Vol. 99 (467)}, 575--583.  

Cressie, N. A. C. (1993). \emph{Statistics for Spatial Data} (revised ed.). Wiley: New York.}

\author{Gel, Y., Raftery, A. E., Gneiting, T., Berrocal, V. J. \email{veronica@stat.washington.edu}.}

\note{Depending on the data, the function might require substantial 
computing time. As a consequence, if the interest is in producing 
probabilistic weather forecasts and generating ensemble members, it is advised to save the output in a file and then use the \code{\link{Variog.fit}} and \code{\link{Field.sim}} functions.  }

\seealso{\code{\link{EmpDir.variog}} for directional empirical variogram averaged over time, and \code{\link{Variog.fit}} for estimation of parameters in a parametric variogram model.}

\examples{
## Loading data
data(slp)
day <- slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <- slp$forecast

## Computing variogram
## No specified cutpoints, no specified maximum distance
## Default number of bins
variogram <- Emp.variog(day=day,obs=obs,forecast=forecast,id=id,
coord1=coord1,coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL)
## Plotting variogram
plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",
ylab="Semi-variance",main="Empirical variogram")

## Computing variogram
## Specified cutpoints, specified maximum distance
## Unspecified number of bins
variogram <- 
Emp.variog(day=day,obs=obs,forecast=forecast,id=id,coord1=coord1,
coord2=coord2,cut.points=seq(0,1000,by=5),max.dist=800,nbins=NULL)
## Plotting variogram
plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",
ylab="Semi-variance",main="Empirical variogram")

\dontshow{## The function is currently defined as
function(day,obs,forecast,id,coord1,coord2,cut.points=NULL,max.dist=NULL,nbins=300){
# default values
if(missing(cut.points))
  cut.points <- NULL
if(missing(max.dist))
  max.dist <- NULL
if(missing(nbins))
  nbins <- NULL
#INPUT CHECK
# Here we check if the input is right.
l.day <- length(day)
l.obs <- length(obs)
l.for <- length(forecast)
l.id <- length(id)
l.coord1 <- length(coord1)
l.coord2 <- length(coord2)
if(sum((c(l.day,l.obs,l.for,l.id,l.coord1,l.coord2)/l.day)==rep(1,6))!=6){
  stop("Error in the dimension of the data")
}
if(sum(is.numeric(obs)==rep("TRUE",l.obs))<l.obs){
  stop("The vector of observations should be a numeric vector")
}
if(sum(ceiling(day)==day)<l.day){
  stop("Day of observation should be an integer")
}
if(sum(is.numeric(coord1)==rep("TRUE",l.coord1)) < l.coord1 | sum(is.numeric(coord2)==rep("TRUE",l.coord2)) < l.coord2){
  stop("Coordinates of the locations should be numeric fields")
}
## here we check the cutpoints vector
l.cuts <- length(cut.points)
if(l.cuts==1){
  stop("Cut points should be a numeric vector")
}

if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))<l.cuts)){
  stop("Cut points should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))==l.cuts)){
  if(sum(order(cut.points)==seq(1:l.cuts)) < l.cuts){
  stop("Cut points should be in increasing order")
  }
  if(sum(cut.points >= 0) < l.cuts){
  stop("Cut points should be non-negative numbers")
  }
  if(length(cut.points)!=length(unique(cut.points))){
  stop("The vector with cut points should not contain repeated entries")
  }
}
   
## check on the max.dist
l.mdist <- length(max.dist)
if(l.mdist > 1){
   print("Max.dist is a numeric field, not a vector")}
if(l.mdist==1){
  if(is.numeric(max.dist)==FALSE){
    stop("Max.dist is a numeric field")
  }
  if(max.dist < 0){
    stop("Max.dist should be a positive number")
  }
}
## check on the number of bins
l.nbins <- length(nbins)
if(l.nbins==0 & l.cuts==0){
  nbins <- 300
}
if(l.nbins==1 & l.cuts >=2){
  nbins=NULL
}
l.nbins <- length(nbins)
if(l.nbins >1){
   stop("Nbins should be an integer: not a vector")
}

if(l.nbins==1){
  if(ceiling(nbins)!=nbins){
     stop("Invalid input: the number of bins should be a positive integer")
  }
  if(ceiling(nbins)==nbins & nbins < 0){
     stop("Invalid input: the number of bins should be a positive integer")
  }
} 
   
  
# ESTIMATION of THE EMPIRICAL VARIOGRAM
# here we order the data in ascending date order
day.o <- order(day)
coord1 <- coord1[day.o]
coord2 <- coord2[day.o]
obs <- obs[day.o]
id <- id[day.o]
forecast <- forecast[day.o]
# the first step in the gop software is to calculate the residuals 
gop.mod <- lm(obs~forecast)
gop.res <- gop.mod$res
gop.var <- var(gop.res)
# the second step is to determine the empirical variograms: if the vector with the cutpoints
# is not specified, we determine the cutpoints by looking at the day with the median average 
# of observations, we calculate the cutpoints so that the number 
# of bins is equal to the one specified and each bin contains approx the same number of pairs. If the vector with the 
# cutpoints is specified, then we just use that vector of cutpoints.

if(length(cut.points)!=0 & length(max.dist)!=0){
  cut.points <- cut.points[cut.points <= max.dist]
}

if(length(cut.points)==0){
# all this part is to determine the day with the median number of observations
  obs.day <- table(day)
  obs.day.vec <- matrix(obs.day,nrow=length(unique(day)),ncol=1)
  obs.median <- round(apply(obs.day.vec,2,median),0)
  diff.median <- abs(obs.day-obs.median)
  unique.day <- unique(day)
  day.index <- min(unique.day[diff.median==min(diff.median)])
# this part is to determine the distance among stations that have been recorded in the 
# day with the median number of observations and to determine the cutpoints.
  obs.day.index <- obs[day==day.index]
  id.day.index <- id[day==day.index]
  coord1.day.index <- coord1[day==day.index]
  coord2.day.index <- coord2[day==day.index]
  dist.day.index <- calc.dist(coord1.day.index,coord2.day.index,id.day.index)
  dist.day.index <- dist.day.index[lower.tri(dist.day.index)]
  dist.day.index <- dist.day.index[dist.day.index > 0]
  ord.dist.day <- sort(dist.day.index)

# this is in case we want to consider only those pairs of stations where the distance 
# is smaller than the maximum distance allowed among locations.
  if(length(max.dist)!= 0){
   ord.dist.day <- ord.dist.day[ord.dist.day < max.dist]
  }
  
  if(length(max.dist)==0){
    max.dist <- quantile(calc.dist(ord.dist.day),.9)
    ord.dist.day <- ord.dist.day[ord.dist.day < max.dist]
  }

  l.dist.day <- length(ord.dist.day)
  l.bins <- floor(l.dist.day/nbins)
  for(i in 1:(nbins-1)){
      cut.points[i] <- ord.dist.day[i*l.bins]
  }
}

# this part is to calculate the empirical variogram
variogram <- avg.variog(day,coord1,coord2,id,gop.res,cut.points)
output <- list(res.var=round(gop.var,3),bin.midpoints=variogram[,1],
              number.pairs=variogram[,2],empir.variog=variogram[,3])
return(output)
  }
}
}
\keyword{file}
