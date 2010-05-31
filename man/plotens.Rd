\name{plotens}
\alias{plotens}
\title{Internal function called by the Field.sim and ProbForecastGOP functions to return a graphical display 
of the members of an ensemble of forecast weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}
\description{
Internal function called by the Field.sim and ProbForecastGOP functions. It returns a graphical display 
of the members of an ensemble of forecast weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}

\usage{
plotens(x,y,grid,number,lims,title,n.pages.4)
}
\arguments{
  \item{x}{numeric vector of length n containing the longitudes of the grid points at which the the members of the ensemble of forecast fields are given.}
  \item{y}{numeric vector of length n containing the latitudes of the grid points at which the the members of the ensemble of forecast fields are given.}
  \item{grid}{numeric array of dimension n by n by number with the values of the members of the ensemble of forecast fields to be displayed (each n by n numeric matrix corresponds to 
  a member of the ensemble of forecast fields to be displayed).}
  \item{number}{number giving the number of members in the ensemble of forecast fields to be displayed.}
  \item{lims}{numeric vector giving the smallest and largest value in the numeric array grid with the values for the ensemble of forecast fields to be displayed.}
  \item{title}{character string giving the title for the plot.}  
  \item{n.pages.4}{number of pages with 4 members to be displayed on screen.}
}
 
\details{
This function is an internal function that is used and called by the Field.sim and the ProbForecastGOP functions to return a graphical display 
of the members of an ensemble of forecast weather field on a region delimited by the lower and upper bound for the longitude and the latitude.
\bold{- Defaults -}
  
None.
}

\value{
The function returns a graphical display 
of the members of an ensemble of forecast weather field on a region delimited by the lower and upper bound for the longitude and the latitude.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





