\name{Field.sim}
\alias{Field.sim}
\title{Simulation of weather random field}
\description{
Simulates and displays realizations of forecast weather fields using a 
geostatistical model. 
}

\usage{
Field.sim(obs, forecast, coord1.grid, coord2.grid, forecast.grid,variog.model="exponential", param.est, n.sim=99, n.displ=4, 
qt.displ=c(10,50,90)) } 

\arguments{
  \item{obs}{numeric vector containing the observed weather quantity.}
  \item{forecast}{numeric vector containing the forecasted weather quantity.}
  \item{coord1.grid}{numeric vector containing the longitudes of the grid points for the forecast.}
  \item{coord2.grid}{numeric vector containing the latitudes of the grid points for the forecast.}
  \item{forecast.grid}{numeric vector containing the forecast grid.}
  \item{variog.model}{character string with the name of the variogram model to be used for the simulations. Implemented models are \kbd{exponential}, \kbd{spherical}, \kbd{gauss}, \kbd{matern}, and \kbd{gencauchy}.}
  \item{param.est}{numeric vector containing values for the parameters to use in the parametric variogram model. 

If the variog.model specified is \kbd{exponential}, \kbd{spherical} or 
\kbd{gauss}, then the parameters required are, in order, the nugget effect 
(a non-negative number), the variance and the range (both positive 
numbers). 

If the variog.model specified is \kbd{gencauchy}, the parameters required 
are, in order, the nugget effect (a non-negative number), the variance, 
the range (both positive numbers), the smoothness parameter \emph{a} (a 
number in \eqn{(0,2]}), and the long-range parameter \emph{b} (a positive 
number).

If the variog.model specified is \kbd{matern} the parameters required are, 
in order, the nugget effect (a non-negative number), the variance, the 
range, and the smoothness parameter \emph{a} (all positive numbers).

For more details on the equation of the variogram models listed above, look below at the section "Details".}
  \item{n.sim}{number of realizations to be simulated.}
  \item{n.displ}{number of realizations to be displayed on screen.}
  \item{qt.displ}{numeric vector containing the quantiles to be displayed.}
}

\details{
The function initially estimates the additive and multiplicative bias by regressing the forecasts on the observed weather quantity. Then it simulates and displays realizations of the weather random field using the specified geostatistical model. Percentiles of the weather random field are calculated using the specified variogram model.

- \bold{Parametric variogram models implemented} - 

Below are the equations of the parametric variogram models that can be implemented.
 
- \kbd{exponential} 

  \deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d}{r}))}

where \eqn{d} is the distance, \eqn{\rho} is the nugget effect, 
\eqn{\sigma^2} is the variance, and \eqn{r} is the range.


- \kbd{spherical}

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (\frac{3}{2}\cdot\frac{d}{r}-\frac{1}{2} \cdot \frac{d^3}{r^3})}

where \eqn{d} is the distance, \eqn{\rho} is the nugget effect, 
\eqn{\sigma^2} is the variance, and \eqn{r} is the range.


- \kbd{gauss}

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d^2}{r^2} ))}

where \eqn{d} is the distance, \eqn{\rho} is the nugget effect, 
\eqn{\sigma^{2}} is the variance, and \eqn{r} is the range.


- \kbd{gencauchy}

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-(1+\frac{d^a}{r^a})^{- \frac{b}{a}})}

where \eqn{d} is the distance, \eqn{\rho} is the nugget effect, 
\eqn{\sigma^{2}} is the variance, \eqn{r} is the range, \eqn{a} is 
the smoothness parameter, and \eqn{b} is the long-range parameter.


- \kbd{matern}

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-\frac{2^{1-a}}{\Gamma(a)}\cdot 
\frac{d^a}{r^a} \cdot K_{a}(\frac{d}{r}))}

where \eqn{d} is the distance, \eqn{\rho} is the nugget effect, 
\eqn{\sigma^{2}} is the variance, \eqn{r} is the range, and \eqn{a} is the 
smoothness parameter.


- Defaults-

By default, \eqn{99} weather random fields are simulated, and \eqn{4} of them are displayed. 

If no vector of percentiles is provided, the \eqn{10}th, \eqn{50}th and 
\eqn{90}th percentiles are determined and displayed. 
}

\value{
The function returns both a graphical and a numerical output. The 
numerical output is a list with the following components:
\item{model}{Character string with the name of the variogram model used.}
\item{nugget}{Value of the nugget effect.}
\item{variance}{Value of the variance.}
\item{range}{Value of the rang.}
\item{additional.par}{Value(s) of the additional parameters required and 
used in the variogram model. This is returned only if the variogram model 
used is \kbd{matern} or \kbd{gencauchy}}
\item{sim.fields}{3-dimensional array where each layer contains the values of the simulated weather field at the gridded locations.}
\item{pct.fields}{3-dimensional array where each layer contains the specified percentile.}

The graphical output consists in plots of the simulated weather fields displayed on multiple pages. 
Before displaying each page, the user will be asked for input.
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: The Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Gneiting, T., Schlather, M. (2004). Stochastic models that separate the 
fractal dimension and the Hurst effect. \emph{SIAM Review} \bold{46}, 
269--282.

Schlather, M. (2001). Simulation and Analysis of Random Fields. 
\emph{R News} \bold{1(2)}, 18--20.

Nychka, D. (2004). \emph{The fields package.} Available at: 
http:lib.stat.cmu.edu/R/CRAN/doc/package/fields.pdf.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data - Some Theory for Kriging}. Springer-Verlag: New York.}

\author{Berrocal, V. J. (\email{veroberrocal@gmail.com}), Raftery, A. E., Gneiting, T., Gel, Y.}

\note{This function depends on the following packages: \pkg{RandomFields} 
and \pkg{fields}.  }

\seealso{\code{\link{plotfields}}, a plotting utility for weather random 
fields, \code{\link[RandomFields]{GaussRF}} in the \pkg{RandomFields} 
package, for simulation of Gaussian random fields.} 

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
coord1.grid <- gridlong$gridded.lon
coord2.grid <- gridlat$gridded.lat
forecast.grid <- forecast.grid$gridded.forecast

## Computing the empirical variogram
variogram <- Emp.variog(day=day,obs=obs,forecast=forecast,id=id,
coord1=coord1,coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL)
## Estimating parameters
param.est <- Variog.fit(emp.variog=variogram,variog.model="exponential",
max.dist.fit=NULL,init.val=NULL,fix.nugget=FALSE)

## Simulating realizations of the weather random field
simul <- 
Field.sim(obs=obs,forecast=forecast,coord1.grid=coord1.grid,
coord2.grid=coord2.grid,forecast.grid=forecast.grid,variog.model="exponential",
param.est=c(param.est$nugget,param.est$variance,param.est$range),n.sim=4,
n.displ=4,qt.displ=c(10,50,90))


\dontshow{
## The function is currently defined as
function(obs,forecast,coord1.grid,coord2.grid,forecast.grid,variog.model="exponential",param.est,n.sim=99,n.displ=4,qt.displ=c(10,50,90)){

## Input check
# set to default values
if(missing(variog.model))
  variog.model <- "exponential"
if(missing(n.sim))
  n.sim <- 99
if(missing(n.displ))
  n.displ <- 4
if(missing(qt.displ))
  qt.displ <- c(10,50,90)
# check on the observations and the forecast
l.obs <- length(obs)
l.for <- length(forecast)
if(l.obs==0){
  stop("obs is a vector")
  }
if(l.for==0){
  stop("forecast is a vector")
  }
if(sum((c(l.obs,l.for)/l.obs)==rep(1,2))!=2){
  stop("Mismatch in dimensions in the data")
  }
if(sum(is.numeric(obs)==rep("TRUE",l.obs))<l.obs){
  stop("obs should be a numeric vector")
  }

if(sum(is.numeric(forecast)==rep("TRUE",l.for))<l.for){
  stop("forecasts should be a numeric vector")
  }

## check on the latitude and forecast grid
l.coord1grid <- length(coord1.grid)
l.coord2grid <- length(coord2.grid)
l.forgrid <- length(forecast.grid)
if(sum((c(l.coord1grid,l.coord2grid,l.forgrid)/l.coord1grid)==rep(1,3))!=3){
  stop("Mismatch in dimensions in the data")
  }
if(sum(is.numeric(coord1.grid)==rep(TRUE,l.coord1grid)) < l.coord1grid){
  stop("Invalid input for coord1.grid")
  }
if(sum(is.numeric(coord2.grid)==rep(TRUE,l.coord2grid)) < l.coord2grid){
  stop("Invalid input for coord2.grid")
  }
if(sum(is.numeric(forecast.grid)==rep(TRUE,l.forgrid)) < l.forgrid){
  stop("Invalid input for forecast.grid")
  }

## check on the variog.model

l.variog <- length(variog.model)
case.var <- 0
if(l.variog==0){
   variog.model <- "exponential"
   case.var <- 1
  }

if(l.variog==1){
  if(is.character(variog.model)=="TRUE"){
    if(variog.model=="exponential"){
        case.var <- 1}
    if(variog.model=="spherical"){
        case.var <- 1}
    if(variog.model=="whittlematern" | variog.model=="matern"){
        case.var <- 1}
    if(variog.model=="gencauchy"){
        case.var <- 1}
    if(variog.model=="gauss"){
        case.var <- 1}
  }
  }

if(case.var==0){
   stop("Incorrect variogram model specification")
  }


## check on the param.est
l.parest <- length(param.est)
case.parest <- 0
if(l.parest < 3){
  stop("Invalid input for param.est")
  }

if(l.parest > 5){
  stop("Invalid input for param.est")
  }


if(l.parest== 3){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<3){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="whittlematern"){
   stop("param.est should be of length 4")
  }
  if(variog.model=="gencauchy"){
   stop("param.est should be of length 5")
  }
  if(sum(param.est >= 0) < 3){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="exponential"){
  case.parest <- 1
  }
  if(variog.model=="spherical"){
  case.parest <- 1
  }
  if(variog.model=="gauss"){
  case.parest <- 1
  }
  }
  

if(l.parest==4){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<4){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="gencauchy"){
   stop("param.est should be of length 5")
  }
  if(variog.model=="exponential"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="gauss"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="spherical"){
   stop("param.est should be of length 3")
  }
  if(sum(param.est >= 0) < 4){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="whittlematern"){
  case.parest <- 1
  }
  }


if(l.parest==5){
  if(sum(is.numeric(param.est)==rep(TRUE,l.parest))<5){
   stop("param.est is a numeric vector")
  }
  if(variog.model=="whittlematern"){
   stop("param.est should be of length 4")
  }
  if(variog.model=="exponential"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="gauss"){
   stop("param.est should be of length 3")
  }
  if(variog.model=="spherical"){
   stop("param.est should be of length 3")
  }
  if(sum(param.est >= 0) < 5){
   stop("The entries of param.est should be non-negative numbers")
  }
  if(variog.model=="gencauchy"){
    if(param.est[4] > 2){
      stop("a is out of the parameter space")
    }
  case.parest <- 1
  }
  }

if(case.parest==0){
  stop("Invalid input for param.est")
  }

## check on the number of simulated random fields

l.sim <- length(n.sim)
if(l.sim==0){
  n.sim <- 99}
l.sim <- length(n.sim)


if(length(n.sim)> 1){
  stop("n.sim is a numeric field: not a vector")
  }

if(length(n.sim)==1){
  if(is.numeric(n.sim)==FALSE){
    stop("The number of simulated fields should be a positive integer")
  }
  if(ceiling(n.sim)!=n.sim){
    stop("The number of simulated fields should be a positive integer")
  }
  if(n.sim <= 0){
    stop("The number of simulated fields should be a positive integer")
  }
  }

## check on the n.displ
l.displ <- length(n.displ)
if(l.displ==0){
  n.displ <- 4}

l.displ <- length(n.displ)
if(l.displ > 1){
  stop("n.displ is a numeric field: not a vector")
  }


if(l.displ==1){
 if(is.numeric(n.displ)!=TRUE){
  stop("n.displ is a numeric field")
 }
 if(ceiling(n.displ)!=n.displ){
  stop("n.displ should be an integer number")
 }
 if(n.displ < 0){
  stop("n.displ should be an integer non-negative number")
 }
 if(n.displ > n.sim){
  stop("n.displ should be less or equal to n.sim")
 }
  }

## check on the qt.displ

l.qt <- length(qt.displ)
if(l.qt==0){
  qt.displ <- c(10,50,90)}
l.qt <- length(qt.displ)

if(l.qt >0){
 if(sum(is.numeric(qt.displ)==rep("TRUE",l.qt))<l.qt){
  stop("qt.displ is a numeric vector")
  }
 if(sum(ceiling(qt.displ)==rep(qt.displ,l.qt))<l.qt){
  stop("The elements of qt.displ should be integer numbers")
  }

 for(i in 1:l.qt){
   if(qt.displ[i] < 0 | qt.displ[i] > 100){
     stop("The elements of qt.displ should be numbers between 0 and 100")
     }
 }
  }


simul <- sim.field(variog.model,param.est,coord1.grid,coord2.grid,n.sim)

# here we add to the 0-mean simulated random fields the bias-corrected gridded forecasts.

gop.mod <- lm(obs~forecast)
gop.coeff <- round(gop.mod$coeff,3)

sim.out <- (gop.mod$coeff[1]+gop.mod$coeff[2]*forecast.grid)+ simul

sim.out.1 <- array(0,c(65,65,n.sim))
for(i in 1:n.sim){
  sim.out.1[,,i] <- engrid(coord1.grid,coord2.grid,sim.out[,i])
  }

# here we determine the percentiles of the random fields
l.qtdispl <- length(qt.displ)
if(l.qtdispl==1 & qt.displ==0){
   quant.out <- 0
   qt.out.1 <- NULL}
  else{
   quant.out <- matrix(0,ncol=l.qtdispl,nrow=l.coord1grid)
   qt.displ <- qt.displ/100
   for(i in 1:l.qtdispl){  
     quant.out[,i] <- (gop.mod$coeff[1]+gop.mod$coeff[2]*forecast.grid) + 
       (qnorm(qt.displ[i],0,1,TRUE,FALSE)*sqrt(param.est[1]+param.est[2]))
   }
   qt.out.1 <- array(0,c(65,65,l.qtdispl))
   for(i in 1:l.qtdispl){
    qt.out.1[,,i] <- engrid(coord1.grid,coord2.grid,quant.out[,i]) 
   }
}

# here we return the output
if(variog.model=="whittlematern" | variog.model=="gencauchy"){
   output <- list(model=variog.model,nugget=param.est[1],variance=param.est[2],
                  range=param.est[3],additional.par=param.est[-seq(1:3)],
                  sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}
if(variog.model=="exponential" | variog.model=="gauss" | variog.model=="spherical"){
   output <- list(model=variog.model,nugget=param.est[1],variance=param.est[2],
                  range=param.est[3],sim.fields=round(sim.out.1,4),pct.fields=round(qt.out.1,4))}
lims <- c(min(sim.out.1[,,1:n.displ],na.rm=TRUE),max(sim.out.1[,,1:n.displ],na.rm=TRUE))

n.displ.4 <- ceiling(n.displ/4)
if(n.displ==1){
     x.lim <- c(min(coord1.grid,na.rm=TRUE),max(coord1.grid,na.rm=TRUE))
     y.lim <- c(min(coord2.grid,na.rm=TRUE),max(coord2.grid,na.rm=TRUE))
     par(ask=TRUE)
     ens.plot(sim.out.1[,,1],lims,x.lim,y.lim,"Ensemble member 1")
}

if(n.displ==2){
     par(mfrow=c(2,2),ask=TRUE)
     plotens(coord1.grid,coord2.grid,sim.out.1[,,1:n.displ],n.displ,lims,"Ensemble member",0)
}

if(n.displ >2){
  if(n.displ.4==1){
     par(mfrow=c(2,2),ask=TRUE)
     plotens(coord1.grid,coord2.grid,sim.out.1[,,1:n.displ],n.displ,lims,"Ensemble member",0)
  }

  if(n.displ.4 >1){
     for(i in 1:n.displ.4){
       n.pages <- i-1 
       n.pages.4 <- 4*n.pages
       if(i!= n.displ.4){
         par(mfrow=c(2,2),ask=TRUE)
         first.col <- (((i-1)*4)+1)
         last.col <- 4*i
         plotens(coord1.grid,coord2.grid,sim.out.1[,,first.col:last.col],4,lims,"Ensemble member",n.pages.4)}
       if(i==n.displ.4){
         par(mfrow=c(2,2),ask=TRUE)
         first.col <- (4*(i-1)+1)
         last.col <- n.displ
         plotens(coord1.grid,coord2.grid,sim.out.1[,,first.col:last.col],((last.col-first.col)+1),lims,"Ensemble member",n.pages.4)}
     }
  }   
}

if(length(qt.out.1)!=0){
   lims <- c(min(qt.out.1[,,1:l.qtdispl],na.rm=TRUE),max(qt.out.1[,,1:l.qtdispl],na.rm=TRUE))
   l.qtdispl.4 <- ceiling(l.qtdispl/4)

   if(l.qtdispl==1){
      x.lim <- c(min(coord1.grid,na.rm=TRUE),max(coord1.grid,na.rm=TRUE))
      y.lim <- c(min(coord2.grid,na.rm=TRUE),max(coord2.grid,na.rm=TRUE))
      title <- paste(qt.displ*100,"-th Percentile")
      par(ask=TRUE)
      ens.plot(qt.out.1[,,1],lims,x.lim,y.lim,title)
   }

   if(l.qtdispl==2){
      par(mfrow=c(2,2),ask=TRUE)
      plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,1:l.qtdispl],l.qtdispl,lims,qt.displ)
   }

   if(l.qtdispl >2){
      if(l.qtdispl.4==1){
        par(mfrow=c(2,2),ask=TRUE)
        plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,1:l.qtdispl],l.qtdispl,lims,qt.displ)
      }

      if(l.qtdispl.4 >1){
        for(i in 1:l.qtdispl.4){
          if(i!= l.qtdispl.4){
            par(mfrow=c(2,2),ask=TRUE)
            first.col <- (((i-1)*4)+1)
            last.col <- 4*i   
            plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,first.col:last.col],4,lims,qt.displ[first.col:last.col])}
     
          if(i==l.qtdispl.4){
            par(mfrow=c(2,2),ask=TRUE)
            first.col <- (4*(i-1)+1)
            last.col <- l.qtdispl  
            plotens.qt(coord1.grid,coord2.grid,qt.out.1[,,first.col:last.col],(last.col-first.col)+1,lims,qt.displ[first.col:last.col])}     
        }   
      }
  }
}
return(output)
}
}
}
\keyword{file}
