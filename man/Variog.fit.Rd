\name{Variog.fit}
\alias{Variog.fit}
\title{Fitting a parametric variogram model to an empirical variogram}
\description{
This function estimates a parametric model fitted to an empirical 
variogram. The estimates are calculated using the weighted least squares 
method.  }

\usage{
Variog.fit(emp.variog, variog.model="exponential", max.dist.fit=NULL, 
init.val=NULL, fix.nugget=FALSE)
}

\arguments{
  \item{emp.variog}{an object of the class \kbd{emp.variog}, output of the function \code{Emp.variog}.}
  \item{variog.model}{character string giving the name of the parametric model to fit to the empirical variogram. Implemented models are \kbd{exponential}, \kbd{spherical}, \kbd{gauss}, \kbd{gencauchy}, \kbd{matern}.}
  \item{max.dist.fit}{number giving the maximum distance considered when fitting the variogram.}
  \item{init.val}{numeric vector giving the initial values for the parameters. The number of initial values to be entered depends on the variogram model specified. 

If the variog.model specified is \kbd{exponential}, \kbd{spherical} or 
\kbd{gauss}, then the parameters required are, in order, the nugget 
effect, the variance and the range. 

If the variog.model specified is \kbd{gencauchy}, the parameters required 
are, in order, the nugget effect, the variance, 
the range, the smoothness parameter a, and the long-range parameter b. 

If the variog.model specified is \kbd{matern} the parameters required are, 
in order, the nugget effect, the variance, the 
range, and the smoothness parameter a. 

For more details on the valid range for the parameters and for the equation of the variogram models listed above, look below at the section "Details".}
  \item{fix.nugget}{logical field indicating whether the nugget should be considered fixed or not. If \kbd{TRUE} the nugget effect will be assumed to be constant, and a value for the fixed nugget effect can be also provided. If the value provided is different from the one entered in the \kbd{init.val} field, then the value of the nugget effect is taken to be the one entered in the \kbd{init.val} field. 

If \kbd{FALSE} the nugget effect will be estimated along with the other parameters in the variogram model.}
}

\details{
The function estimates the parameters in the variogram model by minimizing the weighted least-square loss function
 
  \deqn{LOSS(\theta) = \sum_k n_k [ \left( \frac{\hat{\gamma_k} - \gamma_k(\theta)}{\gamma_k(\theta)} \right)^2 ]}

where \eqn{n_k} is the number of pairs contributing to the variogram computation in the \eqn{k}th bin, \eqn{\hat{\gamma_k}} is the value of the empirical variogram in the \eqn{k}th bin and \eqn{\gamma_k(\theta)} is the value of the estimated parametric variogram model at the midpoint of the \eqn{k}th bin.

- Parametric variogram models -

The parametric model implemented for the variogram are: exponential (\kbd{exponential}), spherical(\kbd{spherical}), gaussian (\kbd{gauss}), generalized Cauchy (\kbd{gencauchy}) and Whittle-Matern (\kbd{matern}).

- \kbd{exponential}:

The equation of the exponential variogram with parameters the nugget effect, \eqn{\rho}, the variance, \eqn{\sigma^{2}}, and the range, \eqn{r}, is given by:

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d}{r}))}

where \eqn{d} is the distance between the two locations and \eqn{\gamma(d)} is the value of the exponential variogram at distance \eqn{d}.

Notice that: the nugget effect, \eqn{\rho}, is a non-negative number, while the variance, \eqn{\sigma^2}, and the range, \eqn{r}, are positive numbers.


- \kbd{spherical}

The equation of the spherical variogram with parameters the nugget effect, \eqn{\rho}, the variance, \eqn{\sigma^{2}}, and the range, \eqn{r}, is given by:

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (\frac{3}{2} 
\cdot\frac{d}{r}-\frac{1}{2} \cdot \frac{d^3}{r^3})}

where \eqn{d} is the distance between the two locations and \eqn{\gamma(d)} is the value of the spherical variogram at distance \eqn{d}. 

Notice that: the nugget effect, \eqn{\rho}, is a non-negative number, while the variance, \eqn{\sigma^2}, and the range, \eqn{r}, are positive numbers.


- \kbd{gauss}

The equation of the gaussian variogram with parameters the nugget effect, \eqn{\rho}, the variance, \eqn{\sigma^{2}}, and the range, \eqn{r}, is given by:

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-exp(- \frac{d^2}{r^2} ))}

where \eqn{d} is the distance between the two locations and \eqn{\gamma(d)} is the value of the gaussian variogram at distance \eqn{d}.

Notice that: the nugget effect, \eqn{\rho}, is a non-negative number, while the variance, \eqn{\sigma^2}, and the range, \eqn{r}, are positive numbers.
 

- \kbd{gencauchy}

The equation of the generalized Cauchy variogram with parameters the nugget effect, \eqn{\rho}, the variance, \eqn{\sigma^{2}}, the range, \eqn{r}, the smoothness parameter, \eqn{a}, and the long-range parameter, \eqn{b}, is given by:

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot 
(1-(1+\frac{d^a}{r^a})^{- \frac{b}{a}})}

where \eqn{d} is the distance between the two locations and \eqn{\gamma(d)} is the value of the generalized Cauchy variogram at distance \eqn{d}.

Notice that: the nugget effect, \eqn{\rho}, is a non-negative number, the 
variance, \eqn{\sigma^2}, and the range, \eqn{r}, are positive numbers, 
the smoothness parameter, \eqn{a}, is a number in \eqn{(0,2]}, and the long-range parameter, \eqn{b}, is a positive number.


- \kbd{matern}

The equation of the Whittle-Matern variogram with parameters the nugget effect, \eqn{\rho}, the variance, \eqn{\sigma^{2}}, the range, \eqn{r}, and the smoothness parameter, \eqn{a}, is given by:

\deqn{\gamma(d) = \rho+\sigma^{2} \cdot (1-(\frac{2^{1-a}}{\Gamma(a)}\cdot 
\frac{d^a}{r^a} \cdot K_{a}(\frac{d}{r}))}

where \eqn{d} is the distance between the two locations, \eqn{\gamma(d)} is the value of the Whittle-Matern variogram at distance \eqn{d}, \eqn{\Gamma} is the gamma function and \eqn{K_a} is the Bessel function of the third kind with parameter \eqn{a}.

Notice that: the nugget effect, \eqn{\rho}, is a non-negative number, the variance, \eqn{\sigma^2}, the range, \eqn{r}, and the smoothness parameter, \eqn{a}, are positive numbers.


- Defaults -

The default value for the maximum distance is \eqn{\frac{1}{2 \cdot \sqrt{2}}} times the maximum distance considered in the empirical variogram. 

Default for the initial values of the parameters is \kbd{NULL}, in which case the initial values for the parameters are determined internally and depend on the empirical variogram values.

By default, the nugget effect is not considered constant. Thus, the 
default value for the \code{fix.nugget} field is \kbd{FALSE}.  
} 

\value{
The function returns a list with components given by: 
  \item{model}{name of the parametric model fitted to the empirical variogram.}
  \item{nugget}{Estimate of the nugget effect.}
  \item{variance}{Estimate of the variance.}
  \item{range}{Estimate of the range.}
  \item{additional.par}{Numeric vector with the estimates of the additional parameters required by the parametric variogram model.}
}

\references{Gel, Y., Raftery, A. E., Gneiting, T. (2004). Calibrated 
probabilistic mesoscale weather field forecasting: the Geostatistical 
Output Perturbation (GOP) method (with discussion). \emph{Journal of the 
American Statistical Association}, \bold{Vol. 99 (467)}, 575--583.

Cressie, N. A. C. (1993) \emph{Statistics for Spatial Data} (revised ed.). Wiley: New York.

Gneiting, T., Schlather, M. (2004). Stochastic models that separate 
fractal dimension and the Hurst effect. \emph{SIAM Review} \bold{46}, 
269--282.

Stein, M. L. (1999). \emph{Interpolation of Spatial Data - Some Theory for Kriging}. Springer-Verlag: New York.}

\author{Gel, Y., Raftery, A. E., Gneiting, T., Berrocal, V. J. \email{veronica@stat.washington.edu}.}

\seealso{\code{\link{Emp.variog}} for computation of the empirical variogram.}

\examples{
## Loading data
library(fields)
library(RandomFields)
data(slp)
day <- slp$date.obs
id <- slp$id.stat
coord1 <- slp$lon.stat
coord2 <- slp$lat.stat
obs <- slp$obs
forecast <- slp$forecast

## Computing the empirical variogram
variogram <- 
Emp.variog(day=day,obs=obs,forecast=forecast,id=id,coord1=coord1,coord2=coord2,
cut.points=NULL,max.dist=NULL,nbins=NULL)
## Estimating parameters
## Without specifying initial values for the parameters
param.variog <- 
Variog.fit(emp.variog=variogram,variog.model="exponential",max.dist.fit=NULL,
init.val=NULL,fix.nugget=FALSE)
## Plotting the empirical variogram with the estimated parametric variogram superimposed
plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",ylab="Semi-variance")
lines(variogram$bin.midpoints,linesmodel(distance=variogram$bin.midpoints,
variog.model="exponential",param=c(param.variog$nugget,
param.variog$variance,param.variog$range)))


## Specifying the initial values for the parameters and keeping the nugget effect fixed
param.variog <- 
Variog.fit(emp.variog=variogram,variog.model="exponential",max.dist.fit=NULL,
init.val=c(0,2,100),fix.nugget=TRUE)
## Plotting the empirical variogram with superimposed the estimated parametric variogram
plot(variogram$bin.midpoints,variogram$empir.variog,xlab="Distance",ylab="Semi-variance")
lines(variogram$bin.midpoints,linesmodel(distance=variogram$bin.midpoints,
variog.model="exponential",param=c(param.variog$nugget,
param.variog$variance,param.variog$range)))


\dontshow{
## The function is currently defined as
function(emp.variog,variog.model="exponential",max.dist.fit=NULL,init.val=NULL,fix.nugget=FALSE){
# INPUT CHECK
### Here there should be the check on whether emp.variog is an object output of the function emp.variog

init.var <- emp.variog$res.var
emp.variog <- matrix(cbind(emp.variog$bin.midpoints,emp.variog$number.pairs,emp.variog$empir.variog),ncol=3)

# Here we do the input check on the rest of the input
# default
if(missing(max.dist.fit))
  max.dist.fit <- NULL
if(missing(fix.nugget))
  fix.nugget <- FALSE

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
   print("Incorrect variogram model specification")
   stop
  }


# check on the max.dist.fit and the initial values

l.max.dist.fit <- length(max.dist.fit)

if(l.max.dist.fit > 1){
  stop("Max.dist.fit should be numeric field: not a vector")
  }

if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==FALSE){
  stop("Max.dist.fit should be a numeric field")
  }

if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==TRUE){
  if(max.dist.fit < 0){
    stop("Max.dist.fit should be positive number")
  }
  if(max.dist.fit > max(emp.variog[,1])){
    stop("Max.dist.fit should be less or equal than max.dist")
  }
  }
 

## Input check on the initial values for the parameters

l.init.val <- length(init.val)

if(l.init.val > 0 & l.init.val <3 ){
  stop("Init.val should be equal to NULL or to a vector of at least length 3")
  }

if(l.init.val ==3 & (sum(is.numeric(init.val)==rep("TRUE",3)) < l.init.val)){
  stop("The initial values should be numeric entries")  
  }


if(l.init.val > 3 & (sum(is.numeric(init.val)==rep("TRUE",3))== l.init.val)){
  if(l.init.val==4 & variog.model!="matern"){
    stop("Incorrect number of initial values")
    }
  if(l.init.val==5 & variog.model!="gencauchy"){
    stop("Incorrect number of initial values")
    }
  if(l.init.val > 5){
    stop("Incorrect number of initial values")
    }
  if(init.val[1] <0){
    stop("The nugget effect cannot be negative")
  }
  if(sum(init.val[2:l.init.val]>0)<(l.init.val-1)){
    stop("Initial values for all the parameters, but the nugget effect, should be positive numbers")
    }
  }


## Input check for the fix.nugget field

l.nug <- length(fix.nugget)
if(l.nug==0){
  fix.nugget <- "FALSE"}
l.nug <- length(fix.nugget)


if(l.nug==1){
  if(fix.nugget!=TRUE & fix.nugget!=FALSE){
   stop("Fix.nugget should be either equal to TRUE or FALSE")
   }
  }

if(l.nug==2){
  if(fix.nugget[1]!=TRUE){
    stop("Invalid input for the fix.nugget field")
  }
  if(is.numeric(fix.nugget[2])==FALSE){
    stop("The second entry of the fix.nugget field should be a numeric field")
  }
  if(fix.nugget[2] <0){
    stop("The second entry of the fix.nugget field should be a non-negative number")
  }
  if(l.init.val >0){
    if(fix.nugget[2]!=init.val[1]){
       fix.nugget[2] <- init.val[1]
    }
  }   
  }

if(l.nug >2){
  stop("Fix.nugget is either a character field or a 2x1 field")
  }

param.est <- round(model.fit(variog.model,emp.variog,max.dist.fit,init.val,init.var,fix.nugget),3)


if(variog.model=="matern"){variog.model <- "whittlematern"}
ifelse((variog.model=="whittlematern" | variog.model=="gencauchy"),
   output <- list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3],additional.par=param.est[-seq(1:3)]),
   output <- list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3]))
   return(output)
  }
}
}
\keyword{file}
