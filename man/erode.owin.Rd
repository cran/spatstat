\name{erode.owin}
\alias{erode.owin}
\title{Erode a Window}
\description{
  Perform morphological erosion of a window 
}
\usage{
 erode.owin(w, r, shrink.frame=TRUE, \dots)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}.}
  \item{r}{positive number: the radius of erosion.}
  \item{shrink.frame}{logical: if \code{TRUE}, erode the bounding
    rectangle as well.}
  \item{\dots}{extra arguments to \code{\link{as.mask}}
    controlling the pixel resolution}
}
\value{
  Another object of class \code{"owin"} representing the
  eroded window.
}
\details{
  The morphological erosion of a set \eqn{W} by a distance \eqn{r > 0}
  is the subset 
  consisting of all points \eqn{x \in W}{x in W} such that the
  distance from \eqn{x} to the boundary of \eqn{W} is greater than
  or equal to \eqn{r}. In other words it is the result of trimming
  a margin of width \eqn{r} off the set \eqn{W}.

  This function computes the erosion of the window \code{w}.
  If \code{w} is not
  a rectangle, it must be approximated by a binary pixel image,
  and the arguments \code{"\dots"} are passed to \code{\link{as.mask}}
  to determine the pixel resolution. There is a sensible default.
  
  If \code{shrink.frame} is false, the resulting window is given the
  same outer, bounding rectangle as the original window \code{w}.
  If \code{shrink.frame} is true, the original bounding rectangle
  is also eroded by the same distance \code{r}.

  To simply compute the area of the eroded window,
  use \code{\link{eroded.areas}}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{eroded.areas}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- erode.owin(w, 0.1) 
  # returns rectangle [0.1, 0.9] x [0.1,0.9]
  \dontrun{
  v <- erode.owin(w, 0.6)
  #fails - erosion is empty
  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
