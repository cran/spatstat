\name{complement.owin}
\alias{complement.owin}
\title{Take Complement of a Window}
\description{
  Take the set complement of a window, within its enclosing rectangle.
}
\usage{
 complement.owin(w)
}
\arguments{
  \item{w}{
    an object of class \code{"owin"} describing
    a window of observation for a point pattern.
  }
}
\value{
  Another object of class \code{"owin"} 
  representing the complement of the window, i.e. the inside
  of the window becomes the outside.
}
\details{
  If \code{w} is a rectangular window this produces an error.
  Otherwise, it yields a window object (of class \code{"owin"},
  see \code{\link{owin.object}}) representing the set complement
  of \code{w} inside its enclosing rectangle. For polygonal windows
  this is achieved by reversing the sign of each boundary polygon,
  while for binary masks it is achieved by negating the pixel values.
}
\seealso{
  \code{\link{owin}},
  \code{\link{owin.object}}
}
\examples{
   # polygonal
   data(demopat)
   w <- demopat$window
   outside <- complement.owin(w)
   # mask
   w <- as.mask(demopat$window)
   outside <- complement.owin(w)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
