\name{setcov}
\alias{setcov}
\title{Set Covariance of a Window}
\description{
  Computes the set covariance function of a window.
}
\usage{
 setcov(W)
}
\arguments{
  \item{W}{A window (object of class \code{"owin"}.}
}
\value{
  A pixel image (an object of class \code{"im"}) representing the
  set covariance function of \code{W}.
}
\details{
  The set covariance function of a region \eqn{W} in the plane
  is the function \eqn{C(v)} defined for each vector \eqn{v}
  as the area of the intersection between \eqn{W} and \eqn{W+v},
  where \eqn{W+v} is the set obtained by shifting (translating)
  \eqn{W} by \eqn{v}.

  We may interpret \eqn{C(v)} as the area of the set of
  all points \eqn{x} in \eqn{W} such that \eqn{x+v} also lies in
  \eqn{W}.
  
  This command computes a discretised approximation to
  the set covariance function of any
  plane region \eqn{W} represented as a window object (of class
  \code{"owin"}, see \code{\link{owin.object}}). The return value is
  a pixel image (object of class \code{"im"}) whose greyscale values
  are values of the set covariance function.

  The set covariance is computed using the Fast Fourier Transform,
  unless \code{W} is a rectangle, when an exact formula is used.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{erode.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  v <- setcov(w)
  plot(v)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
