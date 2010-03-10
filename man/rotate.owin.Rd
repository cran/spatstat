\name{rotate.owin}
\alias{rotate.owin}
\title{Rotate a Window}
\description{
  Rotates a window
}
\usage{
 \method{rotate}{owin}(X, angle=pi/2, \dots)
}
\arguments{
  \item{X}{A window (object of class \code{"owin"}).}
  \item{angle}{Angle of rotation.}
  \item{\dots}{Ignored.}
}
\value{
  Another object of class \code{"owin"} representing the
  rotated window.
}
\details{
  Rotates the window by the specified angle. Angles are measured in
  radians, anticlockwise. The default is to rotate the window 90 degrees
  anticlockwise. The centre of rotation is the origin.

  Rotation of binary image masks is not yet implemented.
}
\seealso{
  \code{\link{owin.object}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- rotate(w, pi/3)
  \dontrun{
  plot(v)
  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
