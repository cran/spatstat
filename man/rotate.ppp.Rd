\name{rotate.ppp}
\alias{rotate.ppp}
\title{Rotate a Point Pattern}
\description{
  Rotates a point pattern
}
\usage{
 \method{rotate}{ppp}(X, angle=pi/2, \dots)
}
\arguments{
  \item{X}{A point pattern (object of class \code{"ppp"}).}
  \item{angle}{Angle of rotation.}
  \item{\dots}{Ignored.}
}
\value{
  Another object of class \code{"ppp"} representing the
  rotated point pattern.
}
\details{
  The points of the pattern, and the window of observation, 
  are rotated about the origin by the angle specified.
  Angles are measured in
  radians, anticlockwise. The default is to rotate the pattern 90 degrees
  anticlockwise. If the points carry marks, these are preserved.

  Rotation of binary image masks is not yet implemented.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{rotate.owin}}
}
\examples{
  data(cells)
  X <- rotate(cells, pi/3)
  \dontrun{
  plot(X)
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
\keyword{math}
