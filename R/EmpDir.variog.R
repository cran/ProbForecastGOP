"EmpDir.variog" <-
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
  stop("Observations should be a numeric vector")
}
if(sum(ceiling(day)==day)<l.day){
  stop("Day of observation should be an integer")
}
if(sum(is.numeric(coord1)==rep("TRUE",l.coord1)) < l.coord1 | sum(is.numeric(coord2)==rep("TRUE",l.coord2)) < l.coord2){
  stop("Coordinates of the locations should be numeric fields")
}
## here we check the tolerance angles

l.tol.angle1 <- length(tol.angle1)

if(l.tol.angle1==0){
  tol.angle1 <- 45}

if(l.tol.angle1 > 1){
  stop("Tol.angle1 should be a number")
}

if(l.tol.angle1==1 & is.numeric(tol.angle1)==FALSE){
  stop("Tol.angle1 should be a number between 0 and 360")
}
if(tol.angle1 <0 | tol.angle1 > 360){
  stop("Tol.angle1 should be a number between 0 and 360")
}

l.tol.angle2 <- length(tol.angle2)

if(l.tol.angle2==0){
  tol.angle2 <- 135}

if(l.tol.angle2 > 1){
  stop("Tol.angle2 should be a number")
}

if(l.tol.angle2==1 & is.numeric(tol.angle2)==FALSE){
  stop("Tol.angle1 should be a number between 0 and 360")
}
if(tol.angle2 <0 | tol.angle2 > 360){
  stop("Tol.angle1 should be a number between 0 and 360")
}

if(tol.angle2 <= tol.angle1){
  stop("Tol.angle2 should be greater than Tol.angle1")
}
## here we check the cutpoints vector
l.cuts <- length(cut.points)
if(l.cuts==1){
  stop("Cutpoints should be a numeric vector")
}
 
if(l.cuts>=2 & (sum(is.numeric(cut.points)==rep("TRUE",l.cuts))<l.cuts)){
  stop("Cutpoints should be a numeric vector")
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
   stop("Max.dist is a numeric field, not a vector")}
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
  stop("Type should be entered")
}
if(type!="E" & type!="N"){
    stop("Type can only be equal to E or to N")
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
# is not specified, we determine the cutpoints by looking at the day with the median average 
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
    max.dist <- quantile(ord.dist.day,.9)
    ord.dist.day <- ord.dist.day[ord.dist.day < max.dist]
  }

  l.dist.day <- length(ord.dist.day)
  l.bins <- ceiling(l.dist.day/nbins)
  for(i in 1:(nbins-1)){
      cut.points[i] <- ord.dist.day[i*l.bins]
  }
}
# this part is to calculate the empirical directional variogram
dir.variogram <- avg.variog.dir(day,coord1,coord2,id,gop.res,cut.points,tol.angle.rad1,tol.angle.rad2,type)
output <- list(bin.midpoints=dir.variogram[,1],
              number.pairs=dir.variogram[,2],dir.variog=dir.variogram[,3])
return(output)
}
