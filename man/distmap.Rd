\name{distmap}
\alias{distmap}
\title{
  Distance Map
}
\description{
  Compute the distance map of an object.
  Generic.
}
\usage{
  distmap(X, \dots)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).
  }
  \item{\dots}{Arguments passed to \code{\link{as.mask}}
    to control pixel resolution.
  }
}
\value{
  A pixel image (object of class \code{"im"}) whose grey scale values
  are the values of the distance map.
}
\details{
  The ``distance map'' of a set of points \eqn{A} is the function
  \eqn{f} whose value \code{f(x)} is defined for any two-dimensional
  location \eqn{x} as the shortest distance from \eqn{x} to \eqn{A}.

  This function computes the distance map of the set \code{X}
  and returns the distance map as a pixel image.
  
  This is generic. Methods are provided for
  point patterns (\code{\link{distmap.ppp}})
  and windows (\code{\link{distmap.owin}}).
}
\seealso{
  \code{\link{distmap.ppp}},
  \code{\link{distmap.owin}}
}
\examples{
  data(cells)
  U <- distmap(cells)
  data(letterR)
  V <- distmap(letterR)
  \dontrun{
  plot(U)
  plot(V)
  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
