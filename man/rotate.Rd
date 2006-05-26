\name{rotate}
\alias{rotate}
\title{Rotate}
\description{
  Applies a rotation to any two-dimensional object,
  such as a point pattern or a window. 
}
\synopsis{
  rotate(X, \dots)
}
\usage{
 rotate(X, angle=pi/2)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{angle}{Angle of rotation.}
}
\value{
  Another object of the same type, representing the
  result of rotating \code{X} through the specified angle.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{rotate.ppp}})
  and windows (\code{\link{rotate.owin}}).

  The angle of rotation is measured in radians, anticlockwise,
  and the centre of rotation is the origin.
}
\seealso{
  \code{\link{rotate.ppp}},
  \code{\link{rotate.owin}}
}
\examples{
  data(cells)
  X <- rotate(cells, pi/3)
  \dontrun{
  plot(X)
  plot(rotate(cells$window, pi/4))
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
