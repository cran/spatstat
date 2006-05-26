\name{dilate.owin}
\alias{dilate.owin}
\title{Dilate a Window}
\description{
  Perform morphological dilation of a window 
}
\usage{
 dilate.owin(w, r, \dots)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}.}
  \item{r}{positive number: the radius of dilation.}
  \item{\dots}{extra arguments to \code{\link{as.mask}}
    controlling the pixel resolution}
}
\value{
  Another object of class \code{"owin"} representing the
  dilated window.
}
\details{
  The morphological dilation of a set \eqn{W} by a distance \eqn{r > 0}
  is the set consisting of all points lying at most \eqn{r} units
  away from \eqn{W}. Effectively, dilation adds a margin of width
  \eqn{r} onto the set \eqn{W}.

  This function computes the dilation of the window \code{w}.
  If \code{w} is not
  a rectangle, it must be approximated by a binary pixel image,
  and the arguments \code{"\dots"} are passed to \code{\link{as.mask}}
  to determine the pixel resolution. There is a sensible default.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- dilate.owin(w, 0.1) 
  # returns rectangle [-0.1, 1.1] x [-0.1, 1.1]
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{math}
