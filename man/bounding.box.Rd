\name{bounding.box}
\alias{bounding.box}
\title{Bounding Box of a Window}
\description{
  Find the smallest rectangle containing a given window.
}
\usage{
 bounding.box(w)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}).}
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"rectangle"} representing a rectangle.
}
\details{
  This function finds the smallest rectangle (with sides parallel to
  the coordinate axes) containing the 
  window represented by \code{w}. 

  Use the function \code{\link{as.rectangle}} to extract the
  bounding window frame, which will be a larger rectangle.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{as.rectangle}}
}
\examples{
  \dontrun{library(spatstat)}
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  r <- bounding.box(w)
  # returns rectangle [1,3] x [2,7]
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
