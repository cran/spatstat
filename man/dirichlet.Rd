\name{dirichlet}
\alias{dirichlet}
\title{Dirichlet Tessellation of Point Pattern}
\description{
  Computes the Dirichlet/Voronoi tessellation of a spatial point pattern.
}
\usage{
dirichlet(X)
}
\arguments{
  \item{X}{Spatial point pattern (object of class \code{"ppp"}).}
}
\details{
  In a spatial point pattern \code{X}, the Dirichlet/Voronoi tile associated
  with a particular point \code{X[i]} is the region of space that is
  closer to \code{X[i]} than to any other point in \code{X}. The
  Dirichlet tiles divide the two-dimensional plane into disjoint
  regions, forming a tessellation.

  This function computes the Dirichlet tessellation (within the original
  window of \code{X}).
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{tess}},
  \code{\link{delaunay}},
  \code{\link{ppp}}
}
\examples{
  X <- runifpoint(42)
  plot(dirichlet(X))
  plot(X, add=TRUE)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
