\name{linesmodel}
\alias{linesmodel}
\title{Computation of parametric variogram model}
\description{
 Computes the value of the parametric variogram model at given distances.
}

\usage{
linesmodel(distance, variog.model="exponential", param)
}

\arguments{
 \item{distance}{numeric vector of distances.}
 \item{variog.model}{character string giving the name of the parametric variogram model. Implemented models are: \kbd{exponential}, \kbd{spherical}, \kbd{gauss}, \kbd{gencauchy} and \kbd{matern}.}
 \item{param}{numeric vector containing the values of the variogram parameters. 
 
If the parametric model specified is \kbd{exponential}, \kbd{spherical} or 
\kbd{gauss}, \code{param} is a vector of length 3 containing, in order: 
the nugget effect (non negative number), the variance and the range (both 
positive numbers).

If the parametric model specified is \kbd{gencauchy}, \code{param} is a 
vector of length 5 whose entries are, in order: the nugget effect (non 
negative number), the variance, the range (both positive numbers), the 
smoothness parameter \kbd{a} (a number in \eqn{(0,2]}), and the long-range 
parameter \kbd{b} (a positive number).

If the parametric model specified is \kbd{matern}, \code{param} is a 
vector of length 4 whose entries are, in order: the nugget effect (a 
non-negative number), the variance, the range, and the smoothness 
parameter \kbd{a} (all three, positive numbers).} }

\details{
 The function calculates the value of the parametric variogram at given distances using the following equations:
 
- If the parametric model is \kbd{exponential}

  \deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d}{r}))}

 where \eqn{\rho} is the nugget effect, \eqn{\sigma^2} is the variance, \eqn{r} is the range, and \eqn{d} is the distance.


- If the parametric model is \kbd{spherical}

  \deqn{\gamma(d) = \rho+\sigma^{2} \cdot 
(\frac{3}{2}\cdot\frac{d}{r}-\frac{1}{2}\cdot \frac{d^3}{r^3})}

 where \eqn{\rho} is the nugget effect, \eqn{\sigma^2} is the variance, \eqn{r} is the range, and \eqn{d} is the distance.


- If the parametric model is \kbd{gauss}

 \deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d^2}{r^2} ))}

 where \eqn{\rho} is the nugget effect, \eqn{\sigma^2} is the variance, \eqn{r} is the range, and \eqn{d} is the distance.


- If the parametric model is \kbd{gencauchy}

 \deqn{\gamma(d) = \rho+\sigma^{2} \cdot 
(1-(1+\frac{d^a}{r^a})^{- \frac{b}{a}})}

 where \eqn{\rho} is the nugget effect, \eqn{\sigma^2} is the variance, 
\eqn{r} is the range, \eqn{d} is the distance, \eqn{a} is the smoothness 
parameter, and \eqn{b} is the long-range parameter.


- If the parametric model is \kbd{matern}

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-(\frac{2^{1-a}}{\Gamma(a)}\cdot 
\frac{d^a}{r^a} \cdot K_{a}(\frac{d}{r}))}

where \eqn{\rho} is the nugget effect, \eqn{\sigma^2} is the variance, 
\eqn{r} is the range, \eqn{d} is the distance, and \eqn{a} is the 
smoothness parameter. }

\value{
The function returns a numeric vector with the values of the parametric variogram model at the bin midpoints.
}

\references{
Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: The Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Cressie, N. A. C. (1993). \emph{Statistics for Spatial Data} (revised ed.). Wiley: New York.

Gneiting, T., Schlather, M. (2004). Stochastic models that separate 
fractal dimension and the Hurst effect. \emph{SIAM Review} \bold{46}, 
269--282.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data - Some Theory for Kriging}. Springer-Verlag: New York.
}

\author{Berrocal, V. J. (\email{veroberrocal@gmail.com}), Raftery, A. E., Gneiting, T., Gel, Y.}

\examples{
## Loading data
data(slp)
day <- slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <- slp$forecast

## Computing empirical variogram
variogram <- Emp.variog(day=day,obs=obs,forecast=forecast,id=id,coord1=coord1,
coord2=coord2,cut.points=NULL,max.dist=NULL,nbins=NULL)

## Estimating variogram parameters
## Without specifying initial values for the parameters
param.variog <- 
Variog.fit(emp.variog=variogram,variog.model="exponential",max.dist.fit=NULL,
init.val=NULL,fix.nugget=FALSE)

## Plotting the empirical variogram with the estimated parametric variogram superimposed
plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",ylab="Semi-variance")
lines(variogram$bin.midpoints,linesmodel(distance=variogram$bin.midpoints,variog.model="exponential",param=c(param.variog$nugget,
param.variog$variance,param.variog$range)))


\dontshow{
## The function is currently defined as
function(distance,variog.model="exponential",param){
# INPUT CHECK
if(missing(distance)){
  stop("Invalid input")
  }

if(missing(variog.model)){
  stop("Invalid input")
  }

if(missing(param)){
  stop("Invalid input")
  }

# Here we check if the vector of distances is well-defined
l.dist <- length(distance)
if(l.dist <1){
  stop("Invalid input")
  }

if(l.dist>=1 & (sum(is.numeric(distance)==rep("TRUE",l.dist)) < l.dist)){
  stop("distance should be a numeric vector")
  }

if(l.dist>=1 & (sum(is.numeric(distance)==rep("TRUE",l.dist))==l.dist)){
  if(sum(distance >=0) < l.dist){
   stop("Entries in distance should be non-negative numbers")
   }
  if(sum(order(distance)==seq(1:l.dist))<l.dist){
   stop("Entries in distance should be in ascending order")
   }
  }

# Here we check the parametric variogram model
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

if(variog.model=="matern"){variog.model <- "whittlematern"}


if(case.var==0){
   stop("Incorrect variogram model specification")
  }


# Here we check the vector with the parameters value
l.par <- length(param)
if(l.par<3){
  stop("Invalid input: the parameter should be a vector of length at least equal to 3")
  }

if(variog.model=="exponential" | variog.model=="spherical" | variog.model=="gauss"){
  if(l.par!=3){
   stop("The parameter vector should have length equal to 3")
  }
  if(l.par==3){
   if(param[1] <0 | param[2] <=0 | param[3] <=0){
    stop("Invalid input for at least one of the parameters")
   }
  }
  }

if(variog.model=="whittlematern"){
  if(l.par!=4){
  stop("The parameter vector should have length equal to 4")
  }
  if(l.par==4){
   if(param[1] <0 | param[2] <=0 | param[3] <=0 | param[4] <=0){
    stop("Invalid input for at least one of the parameters")
   }
  }
  }


if(variog.model=="gencauchy"){
  if(l.par!=5){
  stop("The parameter vector should have length equal to 5")
  }
  if(l.par==5){
   if(param[1] <0 | param[2] <=0 | param[3] <=0 | (param[4] <=0 & param[4]>2) | param[5] <= 0){
    stop("Invalid input for at least one of the parameters")
   }
  }
  }



# This is where the function really starts

if(variog.model=="spherical"){
 fitted.values <- (param[1]+param[2]*(1.5*(distance/param[3])-0.5*(distance/param[3])^3))}
    
if(variog.model=="gauss"){
 fitted.values <- (param[1]+param[2]*(1-exp(-(distance/param[3])^2)))}

if(variog.model=="exponential"){
 fitted.values <- (param[1]+param[2]*(1-exp(-distance/param[3])))}

if(variog.model=="whittlematern" | variog.model=="matern"){
 fitted.values <- 
(param[1]+param[2]-param[2]*(((2^(1-param[4]))/gamma(param[4]))*((distance/param[3])^(param[4]))*besselK(distance/param[3],param[4],expon.scaled=FALSE)))}

if(variog.model=="gencauchy"){
 fitted.values <- (param[1]+param[2]*(1-(1+((distance^2)/param[3])^param[4])^(-param[5]/param[4])))}

return(fitted.values)
  }
}
}
\keyword{file}
