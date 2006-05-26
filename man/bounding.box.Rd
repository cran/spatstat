\name{bounding.box}
\alias{bounding.box}
\title{Bounding Box of a Window}
\description{
  Find the smallest rectangle containing a given window(s)
}
\usage{
 bounding.box(\dots)
}
\arguments{
  \item{\dots}{One or more windows (objects of class \code{"owin"}).}
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"rectangle"} representing a rectangle.
}
\details{
  If the argument is a single window object \code{w},
  this function finds the smallest rectangle (with sides parallel to
  the coordinate axes) containing \code{w}.

  If the argument is a list of windows, then
  this function finds the smallest rectangle that contains
  all of the windows.

  Note that the bounding box of a window \code{w} is usually smaller than
  the window frame (determined by \code{w$xrange} and
  \code{w$yrange} and extracted by \code{\link{as.rectangle}}).
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.rectangle}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- bounding.box(w)
  # returns rectangle [1,3] x [2,7]

  w2 <- unit.square()
  r <- bounding.box(w, w2)
  # returns rectangle [0,3] x [0,7]
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{utilities}
