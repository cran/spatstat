\name{distmap.ppp}
\alias{distmap.ppp}
\title{
  Distance Map of Point Pattern
}
\description{
  Computes the distance from each pixel to the nearest
  point in the given point pattern.
}
\usage{
  distmap.ppp(X, \dots)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{Arguments passed to \code{\link{as.mask}}
    to control pixel resolution.
  }
}
\value{
  A pixel image (object of class \code{"im"}) whose greyscale values
  are the values of the distance map.
  The return value has attributes \code{"index"} and \code{"bdry"}
  which are also pixel images. 
}
\details{
  The ``distance map'' of a point pattern \eqn{X} is the function
  \eqn{f} whose value \code{f(u)} is defined for any two-dimensional
  location \eqn{u} as the shortest distance from \eqn{u} to \eqn{X}.

  This function computes the distance map of the point pattern \code{X}
  and returns the distance map as a pixel image. The greyscale value
  at a pixel \eqn{u} equals the distance from \eqn{u} to the nearest
  point of the pattern \code{X}.

  Additionally, the return value 
  has two attributes, \code{"index"} and \code{"bdry"}, which are
  also pixel images. The grey values in \code{"bdry"} give the
  distance from each pixel to the bounding rectangle of the image.
  The grey values in \code{"index"} are integers identifying which
  point of \code{X} is closest. 
  
  This is a method for the generic function \code{\link{distmap}}.
}
\seealso{
  \code{\link{distmap}},
  \code{\link{distmap.owin}}
}
\examples{
  data(cells)
  U <- distmap(cells)
  \dontrun{
  plot(U)
  plot(attr(U, "bdry"))
  plot(attr(U, "index"))
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
