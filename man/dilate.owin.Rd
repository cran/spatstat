\name{dilate.owin}
\alias{dilate.owin}
\alias{dilation.owin}
\title{Dilate a Window}
\description{
  Perform morphological dilation of a window 
}
\usage{
 dilation.owin(w, r, \dots, polygonal=TRUE, tight=TRUE)
 dilate.owin(w, r, \dots, polygonal=TRUE, tight=TRUE)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}.}
  \item{r}{positive number: the radius of dilation.}
  \item{\dots}{extra arguments to \code{\link{as.mask}}
    controlling the pixel resolution, if the pixel approximation is
    used.
  }
  \item{polygonal}{
    Logical flag indicating whether to compute a polygonal
    approximation to the erosion (\code{polygonal=TRUE}) or
    a pixel grid approximation (\code{polygonal=FALSE}).
  }
  \item{tight}{
    Logical flag indicating whether the bounding frame of the window
    should be taken as the smallest rectangle enclosing the dilated region
    (\code{tight=TRUE}), or should be the
    dilation of the bounding frame of \code{w} (\code{tight=FALSE}).
  }
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

  The functions \code{dilate.owin} and \code{dilation.owin} are
  identical, and compute the dilation of the window \code{w}.

  If \code{polygonal=TRUE} and \code{w} is a rectangle or a polygonal
  window, then a polygonal approximation to the dilation is computed.

  Otherwise, the window \code{w} is first approximated by a binary pixel image,
  and the arguments \code{"\dots"} are passed to \code{\link{as.mask}}
  to determine the pixel resolution. There is a sensible default.
}
\seealso{
  \code{\link{erode.owin}} for the opposite operation.
  
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- dilate.owin(w, 0.1) 
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
