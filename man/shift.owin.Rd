\name{shift.owin}
\alias{shift.owin}
\title{Apply Vector Translation To Window}
\description{
  Applies a vector shift to a window
}
\usage{
 shift.owin(X, vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Window (object of class \code{"owin"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another window (of class \code{"owin"}) representing the
  result of applying the vector shift.
}
\details{
  The window is translated by the vector \code{vec}.
  This is a method for the generic function \code{\link{shift}}.
}
\seealso{
  \code{\link{shift}},
  \code{\link{shift.ppp}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  \dontrun{library(spatstat)}
  W <- owin(c(0,1),c(0,1))
  X <- shift(W, c(2,3))
  \dontrun{
  plot(W)
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
