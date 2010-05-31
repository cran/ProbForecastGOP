\name{engrid}
\alias{engrid}
\title{Internal function used to generate a numeric square matrix that is subsequently used to produce graphical displays of the weather field.}
\description{
Generates a numeric square matrix containing the values of the weather field on a grid to be displayed. Internal function called by the plotfields function.}

\usage{
engrid(coord1, coord2, val)
}
\arguments{
  \item{coord1}{numeric vector containing the longitudes of the meteorological stations.}
  \item{coord2}{vector containing the latitudes of the metereological stations.}
  \item{val}{numeric vector with the values of the weather field to display.}
}
 
\details{
This function is an internal function that is used and called by the plotfields functions to generate a numeric square matrix subsequently used to produce graphical displays of the weather field.

\bold{- Defaults -}
  
None.
}

\value{
The function returns a numeric square matrix with values of the random variables evaluated over a square grid placed on the domain
delimited by the coord1 and coord2 vectors.}

\author{Berrocal, V. J. \email{veroberrocal@gmail.com}, Gel, Y., Raftery, A. E., Gneiting, T.}

\keyword{internal}    





