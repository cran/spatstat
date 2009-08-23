\name{areaLoss}
\alias{areaLoss}
\title{Difference of Disc Areas}
\description{
  Computes the area of that part of a disc
  that is not covered by other discs.
}
\usage{
   areaLoss(X, r, ngrid=spatstat.options("ngrid.disc"), subset=NULL)
}
\arguments{
  \item{X}{
    Locations of the centres of discs.
    A point pattern (object of class \code{"ppp"}).
  }
  \item{r}{
    Disc radius, or vector of disc radii.
  }
  \item{ngrid}{
    Dimensions of the \code{ngrid * ngrid} array of test points
    that will be used to estimate the uncovered area.
  }
  \item{subset}{
    Optional. Index identifying a subset of the points of \code{X}
    for which the area difference should be computed.
  }
}
\value{
  A matrix with one row for each point in \code{X} (or \code{X[subset]})
  and one column for each value in \code{r}. 
}
\details{
  This function computes, for each point \code{X[i]} in \code{X}
  and for each radius \code{r},
  the area of that part of the disc of radius \code{r} centred at the
  location \code{X[i]} that is \emph{not} covered by any of the
  other discs of radius \code{r} centred at the points \code{X[j]}
  for \code{j} not equal to \code{i}.
  This area is important in some calculations related to
  the area-interaction model \code{\link{AreaInter}}.

  The result is a matrix, with one row for each point in \code{X}
  and one column for each entry of \code{r}. 

  The area is calculated by discrete approximation, using an
  \code{ngrid * ngrid} grid of test points. 
}
\seealso{
  \code{\link{AreaInter}},
  \code{\link{areaGain}}
}
\examples{
   data(cells)
   areaLoss(cells, 0.1, 32)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
