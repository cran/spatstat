\name{bdist.pixels}
\alias{bdist.pixels}
\title{Distance to Boundary of Window}
\description{
  Computes the distances 
  from each pixel in a window to the boundary of the window.
}
\usage{
 bdist.pixels(w, \dots, coords=TRUE)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}).}
  \item{\dots}{Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.}
  \item{coords}{Logical: if \code{TRUE}, the result can be
    displayed using \code{persp}, \code{contour} etc.}
}
\value{
  If \code{coords} is false,
  a matrix giving the distances from each pixel in the image raster
  to the boundary of the window. Rows of this matrix correspond to
  the \eqn{y} coordinate and columns to the \eqn{x} coordinate.
  
  If \code{coords} is true, a list with three components
  \code{x,y,z}, where \code{x,y} are vectors of length \eqn{m,n}
  giving the \eqn{x} and \eqn{y} coordinates respectively,
  and \code{z} is an \eqn{m \times n}{m x n} matrix such that
  \code{z[i,j]} is the distance from \code{(x[i],y[j])} to the
  boundary of the window. Rows of this matrix correspond to the
  \eqn{x} coordinate and columns to the \eqn{y} coordinate.
  This result can be plotted with \code{persp}, \code{image}
  or \code{contour}.
}
\details{
  This function computes, for each pixel \eqn{u}
  in the window \code{w}, the shortest distance
  \eqn{d(u, W^c)}{dist(u, W')} from \eqn{u}
  to the boundary of \eqn{W}.

  If the window is not of type \code{"mask"} then it is first
  converted to that type. The arguments \code{"\dots"} are
  passed to \code{\link{as.mask}} to determine the pixel resolution.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{erode.owin}},
  \code{\link{bdist.points}},
}
\examples{
  u <- owin(c(0,1),c(0,1))
  d <- bdist.pixels(u, eps=0.01)
  image(d)
  d <- bdist.pixels(u, eps=0.01, coords=FALSE)
  mean(d >= 0.1)
  # value is approx (1 - 2 * 0.1)^2 = 0.64
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
