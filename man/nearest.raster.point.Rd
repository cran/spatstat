\name{nearest.raster.point}
\alias{nearest.raster.point}
\title{Find Pixel Nearest to a Given Point}
\description{
  Given cartesian coordinates, find the nearest pixel.
}
\usage{
 nearest.raster.point(x,y,w)
}
\arguments{
  \item{x}{Numeric vector of \eqn{x} coordinates of any points}
  \item{y}{Numeric vector of \eqn{y} coordinates of any points}
  \item{w}{A window (an object of class \code{"owin"})
    of type \code{"mask"} representing a binary pixel image.
  }
}
\value{
  List containing two vectors \code{rr} and \code{cc}
  giving row and column positions (in the image matrix).
}
\details{
  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details).
  A window of type \code{"mask"} represents a binary pixel image.

  This function returns a list containing two vectors \code{rr} and
  \code{cc} giving row and column positions (in the image matrix).
  For the location \code{(x[i],y[i])} the nearest
  pixel is at row \code{rr[i]} and column \code{cc[i]} of
  the image.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.mask}}
}
\examples{
  \dontrun{library(spatstat)}
  w <- owin(c(0,1), c(0,1), mask=matrix(TRUE, 100,100))  # 100 x 100 grid
  nearest.raster.point(0.5, 0.3, w)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
