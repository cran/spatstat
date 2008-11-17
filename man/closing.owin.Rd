\name{closing.owin}
\alias{closing.owin}
\title{Morphological Closing of a Window}
\description{
  Perform morphological closing of a window 
}
\usage{
 closing.owin(w, r, \dots)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}.}
  \item{r}{positive number: the radius of the closing.}
  \item{\dots}{extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution}
}
\value{
  Another object of class \code{"owin"} representing the
  closed window.
}
\details{
  The morphological closing (Serra, 1982)
  of a set \eqn{W} by a distance \eqn{r > 0}
  is the set of all points that cannot be
  separated from \eqn{W} by any circle of radius \eqn{r}.
  That is, a point \eqn{x} belongs to the closing \eqn{W*}
  if it is impossible to draw any circle of radius \eqn{r} that
  has \eqn{x} on the inside and \eqn{W} on the outside.
  The closing \eqn{W*} contains the original set \eqn{W}.

  For a small radius \eqn{r}, the closing operation
  has the effect of smoothing out irregularities in the boundary of
  \eqn{W}. For larger radii, the closing operation smooths out
  concave features in the boundary. For very large radii,
  the closed set \eqn{W*} becomes more and more convex.

  This function computes the closing of the window \code{w}
  as a binary pixel mask. If \code{w} is not already a mask, it is first
  converted to a mask by \code{\link{as.mask}}. The arguments
  \code{"\dots"} determine the pixel resolution. There is a sensible
  default.

  The algorithm simply applies \code{\link{dilate.owin}} followed by
  \code{\link{erode.owin}}. 
}
\seealso{
  \code{\link{opening.owin}} for the opposite operation.

  \code{\link{dilate.owin}}, \code{\link{erode.owin}} for the basic
  operations.  
  
  \code{\link{owin}},
  \code{\link{as.owin}} for information about windows.
}
\examples{
  data(letterR)
  v <- closing.owin(letterR, 0.25, dimyx=256)
  plot(v, main="closing.owin")
  plot(letterR, add=TRUE)
}
\references{
  Serra, J. (1982)
  Image analysis and mathematical morphology.
  Academic Press.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
