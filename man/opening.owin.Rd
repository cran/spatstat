\name{opening.owin}
\alias{opening.owin}
\title{Morphological Opening of a Window}
\description{
  Perform morphological opening of a window 
}
\usage{
 opening.owin(w, r, \dots)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}.}
  \item{r}{positive number: the radius of the opening.}
  \item{\dots}{extra arguments passed to \code{\link{as.mask}}
    controlling the pixel resolution}
}
\value{
  Another object of class \code{"owin"} representing the
  opened window.
}
\details{
  The morphological opening (Serra, 1982)
  of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset of points in \eqn{W} that can be 
  separated from the boundary of \eqn{W} by a circle of radius \eqn{r}.
  That is, a point \eqn{x} belongs to the opening 
  if it is possible to draw a circle of radius \eqn{r} (not necessarily
  centred on \eqn{x}) that has \eqn{x} on the inside
  and the boundary of \eqn{W} on the outside.
  The opened set is a subset of \code{W}.

  For a small radius \eqn{r}, the opening operation
  has the effect of smoothing out irregularities in the boundary of
  \eqn{W}. For larger radii, the opening operation removes promontories
  in the boundary. For very large radii, the opened set is empty.

  This function computes the opening of the window \code{w}
  as a binary pixel mask. If \code{w} is not already a mask, it is first
  converted to a mask by \code{\link{as.mask}}. The arguments
  \code{"\dots"} determine the pixel resolution. There is a sensible
  default.

  The algorithm simply applies \code{\link{erode.owin}} followed by
  \code{\link{dilate.owin}}. 
}
\seealso{
  \code{\link{closing.owin}} for the opposite operation.

  \code{\link{dilate.owin}}, \code{\link{erode.owin}} for the basic
  operations.  
  
  \code{\link{owin}},
  \code{\link{as.owin}} for information about windows.
}
\examples{
  data(letterR)
  v <- opening.owin(letterR, 0.3, dimyx=256)
  plot(v, main="opening.owin")
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
