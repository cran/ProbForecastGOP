\name{ens.plot}
\alias{ens.plot}
\title{Internal function used to return a graphical display 
of a weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}
\description{
Internal function called by the ProbForecastGOP function. It returns a graphical display 
of a weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}

\usage{
ens.plot(grid,lims,x.lim,y.lim,title)
}
\arguments{
  \item{grid}{numeric matrix of dimension n by n with the values of the ensemble member to be displayed.}
  \item{lims}{numeric vector giving the smallest and the largest value of the ensemble member to be displayed.}
  \item{x.lim}{numeric vector giving the smallest and the largest longitude to be displayed.}
  \item{y.lim}{numeric vector giving the smallest and the largest latitude to be displayed.}
  \item{title}{character string with the title for the plot.}  
}
 
\details{
This function is an internal function that is used and called by the ProbForecastGOP function to returns a graphical display 
of a weather field on a region delimited by the lower and upper bound for the longitude and the latitude.
\bold{- Defaults -}
  
None.
}

\value{
The function returns a graphical display of a weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





