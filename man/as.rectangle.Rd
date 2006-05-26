\name{as.rectangle}
\alias{as.rectangle}
\title{Window Frame}
\description{
  Get the window frame of an (irregular) window
}
\usage{
 as.rectangle(\dots)
}
\arguments{
  \item{\dots}{A window. Either an object of class \code{"owin"},
    or other data determining such a window. The argument(s)
    are passed to \code{\link{as.owin}}.}
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"rectangle"} representing a rectangle.
}
\details{
  This function just extracts the outer bounding rectangle
  of the window \code{w} as given by its elements \code{xrange,yrange}.

  Use the function \code{\link{bounding.box}} to compute the smallest
  bounding rectangle.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{bounding.box}}
}
\examples{
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- as.rectangle(w)
  # returns a 10 x 10 rectangle
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{manip}
