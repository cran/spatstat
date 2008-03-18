\name{areadiff}
\alias{areadiff}
\title{Difference of Disc Areas}
\description{
  Computes the area of that part of a disc
  that is not covered by other discs.
}
\usage{
   areadiff(u, X, r, ngrid=256)
}
\arguments{
  \item{u}{
    Coordinates of the centre of the disc.
    A vector of length 2.
  }
  \item{X}{
    Locations of the centres of other discs.
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Disc radius.
  }
  \item{ngrid}{
    Dimensions of the \code{ngrid * ngrid} array of test points
    that will be used to estimate the uncovered area.
  }
}
\value{
  A single number (the uncovered area)
}
\details{
  This function computes the area of that part of
  the disc of radius \code{r} centred at the location \code{u}
  that is \emph{not} covered by any of the discs of radius \code{r}
  centred at the points of the pattern \code{X}. 
}
\seealso{
  \code{\link{AreaInter}}
}
\examples{
   data(cells)
   u <- c(0.5,0.5)
   areadiff(u, cells, 0.1)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
