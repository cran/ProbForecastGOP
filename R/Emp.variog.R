"Emp.variog" <-
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
   stop("Max.dist is a numeric field, not a vector")
}
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
    max.dist <- quantile(ord.dist.day,.9)
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
