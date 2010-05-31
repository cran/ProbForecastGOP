\name{model.fit}
\alias{model.fit}
\title{Computation of parametric variogram model}
\description{
 Internal function used and called by the Variog.fit and the ProbForecastGOP functions to fit a parametric variogram model to an empirical variogram.
}

\usage{
model.fit(par.model,emp.variog,max.dist.fit,init.val,init.var,fix.nugget)
}

\arguments{
  \item{par.model}{character string giving the name of the parametric 
model to fit to the empirical variogram. Implemented models are 
\kbd{exponential}, \kbd{spherical}, \kbd{gauss}, \kbd{gencauchy}, and 
\kbd{matern}.}
  \item{emp.variog}{an object of the class \kbd{emp.variog}, output of the function \code{Emp.variog} or of the function 
\code{avg.variog}.}
  \item{max.dist.fit}{number giving the maximum distance considered when fitting the variogram.}
  \item{init.val}{numeric vector giving the initial values for the parameters. The number of initial values to be entered depends on the variogram model specified. }
  \item{init.var}{number giving an initial estimate for the variance parameter of the variogram.}
  \item{fix.nugget}{logical field indicating whether the nugget should be considered fixed or not. If \kbd{TRUE} the nugget effect will be assumed to be constant, and a value 
  for the fixed nugget effect can be also provided. If the value provided is different from the one entered in the \kbd{init.val} field, then the value of the nugget effect 
  is taken to be the one entered in the \kbd{init.val} field. }
}

\details{
The function estimates the parameters in the variogram model by minimizing the weighted least-square loss function.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a numeric vector with the estimates of the parameters of the parametric variogram model fitted to the empirical variogram.
}


\author{Berrocal, V. J. (\email{veroberrocal@gmail.com}), Raftery, A. E., Gneiting, T., Gel, Y.}
\keyword{internal}
