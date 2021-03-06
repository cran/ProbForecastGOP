\name{EmpDir.variog}
\alias{EmpDir.variog}
\title{Directional empirical variogram of forecast errors averaged over time}
\description{
Calculates directional empirical variogram of forecast errors, averaged over time. 
}

\usage{EmpDir.variog(day, obs, forecast, id, coord1, coord2, tol.angle1=45, tol.angle2=135, cut.points=NULL, max.dist=NULL, nbins=300, type) }

\arguments{
  \item{day}{numeric vector containing the day of observation.}
  \item{obs}{numeric vector containing the observed weather quantity.}
  \item{forecast}{numeric vector containing the forecasted weather quantity.}
  \item{id}{vector with the id of the metereological stations.}
  \item{coord1}{vector containing the longitudes of the metereological stations.}
  \item{coord2}{vector containing the latitudes of the metereological stations.}
  \item{tol.angle1}{number giving a lower bound for the tolerance angle (measured in degrees).}
  \item{tol.angle2}{number giving an upper bound for the tolerance angle (measured in degrees).}
  \item{cut.points}{numeric vector containing the cutpoints used for variogram binning.}
  \item{max.dist}{a numerical value giving the upper bound for the distance considered in the variogram computation}.
  \item{nbins}{a numerical value giving the number of bins to use for variogram binning. If both \code{cut.points} and \code{nbins} are entered, the entry for \code{nbins} 
will be ignored and the vector with the cutpoints will instead be used for variogram binning.}
  \item{type}{character string indicating the direction to use for variogram computations. Possible values are either 'E' (for East-West) or 'N' (for North-South).}
}

\details{ 
The function includes bias-correction; it regresses the forecasts on the observed weather quantity and computes the residuals. The directional empirical variogram 
of the residuals is then calculated by determining, for each day, the "directional" distance among all pairs of stations that have been observed in the same day and by 
calculating for each day the sum of all the squared differences in the residuals within each bin. These sums are then averaged over time, with weights for each bin given by 
the sum over time of the number of pairs of stations within the bin.

The formula used is:
 \deqn{ \gamma(h) = \sum_d \frac{1}{2N_{(h,d)}} (\sum_i (Y(x_{i}+h,d)-Y(x_{i},d))^2)}

where \eqn{\gamma(h)} is the empirical variogram at distance \eqn{h}, \eqn{N_{(h,d)}} is the number of pairs of stations that have been recorded at day \eqn{d} and whose 
distance is equal to \eqn{h}, and \eqn{Y(x_{i}+h,d)} and \eqn{Y(x_{i},d)} are, respectively, the values of the residuals on day \eqn{d} at stations \eqn{x_{i}+h} and 
\eqn{x_{i}}. Variogram binning is ignored in this formula. 

The "directional" distance between two locations is defined to be equal to the distance between the two locations if the angle between the two locations is within the 
allowed range, while it is set equal to infinity if the angle between the two locations is outside the allowed range.

\bold{- Defaults -}

By default, \code{tol.angle1} and \code{tol.angle2} are set to 45 and 135 degrees, respectively. If the vector with the cutpoints is not specified, the cutpoints are 
determined so that there are \code{nbins} number of bins with approximately the same number of pairs per bin.

If both the vector with the cutpoints and the number of bins, \code{nbins}, are not provided, the function by default determines the cutpoints so that there are a total of 300 bins with approximately the same number of pairs per bin. If both the vector with the cutpoints and the number of bins are provided, the entry for the number of bins is ignored and the vector with the cutpoints is used for variogram binning. 

The default value for the maximum distance considered in the variogram computation is the \eqn{90}-th percentile of the distances between the stations.  
}

\value{
The function returns a list with components given by:
  \item{bin.midpoints}{Numeric vector with midpoints of the bins used in the directional empirical variogram computation.}
  \item{number.pairs}{Numeric vector with the number of pairs per bin.}
  \item{dir.variog}{Numeric vector with the directional empirical variogram values.}
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: The Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Gel, Y., Raftery, A. E., Gneiting, T., Berrocal, V. J. (2004). 
Rejoinder. \emph{Journal of the American Statistical Association}, 
\bold{Vol. 99 (467)}, 588--590.

Cressie, N. A. C. (1993). \emph{Statistics for Spatial Data} (revised ed.). Wiley: New York.
}

\author{Berrocal, V. J. (\email{veroberrocal@gmail.com}), Raftery, A. E., Gneiting, T., Gel, Y.}

\note{The function might require some time to return an output.}

\seealso{\code{\link{Emp.variog}} for empirical variogram averaged over time, \code{\link{avg.variog}} and \code{\link{avg.variog.dir}} for, respectively, empirical and 
directional empirical variogram of a random variable averaged over time, and \code{\link{Variog.fit}} for estimation of parameters in a parametric variogram model. }

\examples{
## Loading data
data(slp)
day <- slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <- slp$forecast

## Computing directional variogram
## No specified cutpoints, no specified maximum distance
## No specified tolerance angles and default number of bins
dir.variog <- EmpDir.variog(day,obs,forecast,id,coord1,coord2,
tol.angle1=NULL,tol.angle2=NULL,cut.points=NULL,max.dist=NULL,
nbins=NULL,type='E')
## Plotting directional variogram
plot(dir.variog$bin.midpoints,dir.variog$dir.variog,xlab="Distance",
ylab="Semi-variance",main="Empirical Directional variogram")


## Computing directional variogram
## Specified cutpoints, specified maximum distance
## Specified tolerance angles and unspecified number of bins
dir.variog <- 
EmpDir.variog(day,obs,forecast,id,coord1,coord2,tol.angle1=30,
tol.angle2=150,cut.points=seq(0,1000,by=5),max.dist=800,nbins=NULL,
type='N')
## Plotting directional variogram
plot(dir.variog$bin.midpoints,dir.variog$dir.variog,xlab="Distance",ylab="Semi-variance",main="Empirical Directional variogram")


\dontshow{## The function is currently defined as
function(day,obs,forecast,id,coord1,coord2,tol.angle1=45,tol.angle2=135,cut.points=NULL,max.dist=NULL,nbins=300,type){

# default values
if(missing(cut.points))
  cut.points <- NULL
if(missing(max.dist))
  max.dist <- NULL
if(missing(nbins))
  nbins <- NULL
if(missing(type)){
  stop("Invalid input for type: a type should be provided")
  }
if(missing(tol.angle1))
  tol.angle1 <- 45
if(missing(tol.angle2))
  tol.angle2 <- 135

#INPUT CHECK
# Here we check if the input is right.
l.day <- length(day)
l.obs <- length(obs)
l.for <- length(forecast)
l.id <- length(id)
l.coord1 <- length(coord1)
l.coord2 <- length(coord2)

if(sum((c(l.day,l.obs,l.for,l.id,l.coord1,l.coord2)/l.day)==rep(1,6))!=6){
  stop("Mismatch in dimensions in the data")
  }
if(sum(is.numeric(obs)==rep("TRUE",l.obs))<l.obs){
  stop("obs should be a numeric vector")
}
if(sum(ceiling(day)==day)<l.day){
  stop("day should be a vector containing integers")
}
if(sum(is.numeric(coord1)==rep("TRUE",l.coord1)) < l.coord1 | sum(is.numeric(coord2)==rep("TRUE",l.coord2)) < l.coord2){
  stop("coord1 and coord2 should be numeric vectors")
}
## here we check the tolerance angles

l.tol.angle1 <- length(tol.angle1)

if(l.tol.angle1==0){
  tol.angle1 <- 45}

if(l.tol.angle1 > 1){
  stop("tol.angle1 should be a number")
}
 
if(l.tol.angle1==1 & is.numeric(tol.angle1)==FALSE){
  stop("tol.angle1 should be a number between 0 and 360")
}

if(tol.angle1 <0 | tol.angle1 > 360){
  stop("tol.angle1 should be a number between 0 and 360")
}

l.tol.angle2 <- length(tol.angle2)

if(l.tol.angle2==0){
  tol.angle2 <- 135}

if(l.tol.angle2 >1){
  stop("tol.angle2 should be a number")
}

if(l.tol.angle2==1 & is.numeric(tol.angle2)==FALSE){
  stop("tol.angle1 should be a number between 0 and 360")
}
if(tol.angle2 <0 | tol.angle2 > 360){
  stop("tol.angle1 should be a number between 0 and 360")
}

if(tol.angle2 <= tol.angle1){
  stop("tol.angle2 should be greater than tol.angle1")
}
## here we check the cutpoints vector
l.cuts <- length(cut.points)
if(l.cuts==1){
  stop("cut.points should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))<l.cuts)){
  stop("cut.points should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))==l.cuts)){
  if(sum(order(cut.points)==seq(1:l.cuts))<l.cuts){
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
   stop("max.dist is a numeric field, not a vector")}
if(l.mdist==1){
  if(is.numeric(max.dist)==FALSE){
    stop("max.dist is a numeric field")
  }
  if(max.dist < 0){
    stop("max.dist should be a positive number")
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
   stop("nbins should be an integer: not a vector")
}
if(l.nbins==1){
  if(floor(nbins)!=nbins){
     stop("Invalid input: the number of bins should be a positive integer")
  }
  if(ceiling(nbins)==nbins & nbins < 0){
     stop("Invalid input: the number of bins should be a positive integer")
  }
} 
   
  
## check on the type ##
l.type <- length(type)
if(l.type==0){
  stop("type should be entered")
}
if(type!="E" & type!="N"){
    stop("type can only be equal to E or to N")
}
# ESTIMATION of THE EMPIRICAL VARIOGRAM
# here we order the data in ascending date order
day.o <- order(day)
coord1 <- coord1[day.o]
coord2 <- coord2[day.o]
obs <- obs[day.o]
id <- id[day.o]
forecast <- forecast[day.o]
# here we calculate the residuals 
gop.mod <- lm(obs~forecast)
gop.res <- gop.mod$res
gop.var <- var(gop.res)
tol.angle.rad1 <- tol.angle1/57.29577951
tol.angle.rad2 <- tol.angle2/57.29577951

# if the vector with the cutpoints
# is not specified, we determine the cutpoints by looking at the day with the median number 
# of observations, we calculate the cutpoints so that the number 
# of bins is equal to the one specified and each bin contains approx the same number of pairs. If the vector with the 
# if the vector of cutpoints is specified, then we just use that vector of cutpoints.

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

# this part is to determine the directional distance among stations that have been recorded in the 
# day with the median number of observations and to determine the cutpoints.
  obs.day.index <- obs[day==day.index]
  id.day.index <- id[day==day.index]
  coord1.day.index <- coord1[day==day.index]
  coord2.day.index <- coord2[day==day.index]
  dist.day.index <- calc.dist.dir(coord1.day.index,coord2.day.index,id.day.index,tol.angle.rad1,tol.angle.rad2,type)
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

# this part is to calculate the empirical directional variogram
unique.day <- unique(day)   
l.day <- length(unique.day)
n.stations <- length(unique(id))
n.cuts <- length(cut.points)
W <- rep(0,(n.cuts-1))
W.obs <- rep(0,(n.cuts-1))
      
for(i in 1:l.day){
      distance.day <- NULL
      difference.day <- NULL
      new.dist.day <- NULL
      s.new.dist.day <- NULL
      o.new.dist.day <- NULL
      new.diff.day <- NULL
      new.diff.day <- NULL
      s.new.diff.day <- NULL
      gop.res.day <- NULL
      id.day <- NULL
      coord1.day <- NULL
      coord2.day <- NULL
      gop.res.day <- gop.res[day==unique.day[i]]
      id.day <- id[day==unique.day[i]]
      coord1.day <- coord1[day==unique.day[i]]
      coord2.day <- coord2[day==unique.day[i]]
      distance.day <-
calc.dist.dir(coord1.day,coord2.day,id.day,tol.angle.rad1,tol.angle.rad2,type)
       new.dist.day <- distance.day[lower.tri(distance.day)]
       s.new.dist.day <- sort(new.dist.day)
       o.new.dist.day <- order(new.dist.day)
       difference.day  <- calc.difference(gop.res.day)
       new.diff.day <- difference.day[lower.tri(difference.day)]
       s.new.diff.day <- new.diff.day[o.new.dist.day]
       W.new <- rep(0,(n.cuts-1))
       W.obs.new <- rep(0,(n.cuts-1))
      
       for(s in 1:(n.cuts-1)){
           low.bound <- cut.points[s]
           upp.bound <- cut.points[s+1]
           v.dist <- NULL
           v.dist1 <- NULL
           v.diff1 <- NULL
           index.v.dist <- NULL
           index.v.dist1 <- NULL
       
           if(s < (n.cuts-1)){
                v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & s.new.dist.day < upp.bound]
                index.v.dist <- seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day < upp.bound]
                v.dist1 <- v.dist[v.dist!=0]
                index.v.dist1 <- index.v.dist[v.dist!=0]
           
                if(length(v.dist1) >=1){
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                if(length(v.dist1) ==0){
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
           if(s==(n.cuts-1)){
                 v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 index.v.dist <- seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 v.dist1 <- v.dist[v.dist!=0]
                 index.v.dist1 <- index.v.dist[v.dist!=0]  
                   
                 if(length(v.dist1) >=1){
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                 if(length(v.dist1) ==0){
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
           if(s==(n.cuts-1)){
                 v.dist <- s.new.dist.day[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 index.v.dist <- seq(1:length(s.new.dist.day))[s.new.dist.day >= low.bound & s.new.dist.day <= upp.bound]
                 v.dist1 <- v.dist[v.dist!=0]
                 index.v.dist1 <- index.v.dist[v.dist!=0]
    
                 if(length(v.dist1) >=1){
                   v.diff1 <- s.new.diff.day[index.v.dist1]
                   W.new[s] <- sum(v.diff1)
                   W.obs.new[s] <- length(v.dist1)}
                 if(length(v.dist1) ==0){
                   W.new[s] <- 0
                   W.obs.new[s] <- 0}
           }
       }
   W <- W+W.new
   W.obs <- W.obs + W.obs.new
   }
   avg.variog <- round(W/(2*W.obs),2)
   n.h <- W.obs
   x.vec <- NULL
   for(i in 1:(n.cuts-1)){
     x.vec[i] <- (cut.points[i]+cut.points[i+1])/2
   }
   fin.avg.variog <- c(0,avg.variog)
   fin.x.vec <- c(0,x.vec)
   fin.n.h <- c(n.stations,n.h)
   B <- list(bin.midpoints=fin.x.vec,number.pairs=fin.n.h,dir.variog=fin.avg.variog)
   return(B)
}
}
}
\keyword{file}
