\name{expand.owin}
\alias{expand.owin}
\title{Expand Window By Factor}
\description{
  Expands the bounding box of a window by a given factor.
}
\usage{
 expand.owin(W, f=1)
}
\arguments{
  \item{W}{A window.}
  \item{f}{Scalar factor for expansion of the \emph{area} of the window.}
}
\value{
  A rectangular window.
}
\details{
  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  or can be given in any format acceptable to \code{\link{as.owin}()}.

  The bounding box of the
  window \code{W} is expanded by a scale factor equal to the square root
  of \code{f} in both the \eqn{x} and \eqn{y} directions,
  so that the area increases by a factor \code{f}. The centre of this
  expanded rectangle is the same as the centre of the original bounding
  box. The expanded rectangle is returned.

  This quirky little function exists mainly for reference,
  as it is the rule for calculating window expansion in our implementation of
  the Metropolis-Hastings algorithm \code{\link{rmh.default}}.

  For general transformations of the scale, location and orientation
  of a window, see \code{\link{shift}}, \code{\link{affine}}
  and \code{\link{rotate}}.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{shift}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  w <- square(1)
  expand.owin(w, 9)
  # returns the square [-1,2] x [-1,2]
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
