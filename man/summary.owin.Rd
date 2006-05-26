\name{summary.owin}
\alias{summary.owin}
\title{Summary of a Spatial Window}
\description{
  Prints a useful description of a window object.
}
\synopsis{
  \method{summary}{owin}(object, \dots)
}
\usage{
 \method{summary}{owin}(object)
}
\arguments{
  \item{object}{Window (object of class \code{"owin"}).}
}
\details{
  A useful description of the window \code{object} is printed.

  This is a method for the generic function \code{\link{summary}}.
}
\seealso{
  \code{\link{summary}},
  \code{\link{summary.ppp}},
  \code{\link{print.owin}}
}
\examples{
  summary(owin())  # the unit square

  data(demopat)
  W <- demopat$window  # weird polygonal window
  summary(W)           # describes it

  summary(as.mask(W))  # demonstrates current pixel resolution
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{methods}

