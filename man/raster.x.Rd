\name{raster.x}
\alias{raster.x}
\alias{raster.y}
\title{Cartesian Coordinates for a Pixel Raster}
\description{
  Returns a matrix giving the \eqn{x} (or \eqn{y}) coordinates
  of each pixel in a binary pixel image.
}
\synopsis{
  raster.x(w)
}
\usage{
 raster.x(w)
 raster.y(w)
}
\arguments{
  \item{w}{A window (an object of class \code{"owin"})
    of type \code{"mask"} representing a binary pixel image.
  }
}
\value{
  A matrix of the same dimensions as the pixel grid in \code{w},
  and giving the value of the \eqn{x} (or \eqn{y}) coordinate
  of each pixel in the raster.
}
\details{
  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details).
  A window of type \code{"mask"} represents a binary pixel image.
  This function returns a matrix of the same dimensions as the
  binary pixel image itself, with entries giving the \eqn{x} coordinate
  (for \code{raster.x}) or \eqn{y} coordinate (for \code{raster.y})
  of each pixel in the image.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.mask}}
}
\examples{
  \dontrun{library(spatstat)}
  u <- owin(c(-1,1),c(-1,1)) # square of side 2
  w <- as.mask(u, eps=0.01) # 200 x 200 grid
  X <- raster.x(w)
  Y <- raster.y(w)
  disc <- owin(c(-1,1), c(-1,1), mask=(X^2 + Y^2 <= 1))
  \dontrun{plot(disc)}
  # approximation to the unit disc
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
