\name{shift}
\alias{shift}
\title{Apply Vector Translation}
\description{
  Applies a vector shift of the plane 
  to a geometrical object,
  such as a point pattern or a window. 
}
\synopsis{
  shift(X, \dots)
}
\usage{
 shift(X, vec=c(0,0))
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{vec}{Vector of length 2 representing a translation.}
}
\value{
  Another object of the same type, representing the
  result of applying the shift.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{shift.ppp}})
  and windows (\code{\link{shift.owin}}).

  The object is translated by the vector \code{vec}.
}
\seealso{
  \code{\link{shift.ppp}},
  \code{\link{shift.owin}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  data(cells)
  X <- shift(cells, c(2,3))
  \dontrun{
  plot(X)
  # no discernible difference except coordinates are different
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
