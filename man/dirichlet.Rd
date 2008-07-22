\name{dirichlet}
\alias{dirichlet}
\title{Dirichlet Tessellation of Point Pattern}
\description{
  Computes the Dirichlet tessellation of a spatial point pattern.
}
\usage{
dirichlet(x)
}
\arguments{
  \item{x}{Spatial point pattern (object of class \code{"ppp"}).}
}
\details{
  In a spatial point pattern \code{x}, the Dirichlet tile associated
  with a particular point \code{x[i]} is the region of space that is
  closer to \code{x[i]} than to any other point in \code{x}. The
  Dirichlet tiles divide the two-dimensional plane into disjoint
  regions, forming a tessellation.

  This function computes the Dirichlet tessellation (within the original
  window of \code{x}).
}
\value{
  A tessellation (object of class \code{"tess"}).
}
\seealso{
  \code{\link{tess}},
  \code{\link{ppp}}
}
\examples{
  X <- runifpoint(42)
  plot(dirichlet(X))
  plot(X, add=TRUE)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
