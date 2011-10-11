\name{imcov}
\alias{imcov}
\title{Spatial Covariance of a Pixel Image}
\description{
  Computes the unnormalised spatial covariance function of a pixel image.
}
\usage{
 imcov(X, Y=X)
}
\arguments{
  \item{X}{
    A pixel image (object of class \code{"im"}.
  }
  \item{Y}{
    Optional. Another pixel image.
  }
}
\value{
  A pixel image (an object of class \code{"im"}) representing the
  spatial covariance function of \code{X},
  or the cross-covariance of \code{X} and \code{Y}.
}
\details{
  The (uncentred, unnormalised)
  \emph{spatial covariance function} of a pixel image \eqn{X} in the plane
  is the function \eqn{C(v)} defined for each vector \eqn{v}
  as the integral 
  \deqn{
    C(u) = \int X(u)X(u+v)\, {\rm d}u
  }{
    C(u) = integral of X(u) * X(u+v) du
  }
  over all spatial locations \eqn{u}, where \eqn{X(u)} denotes the
  pixel value at location \eqn{u}.
  
  This command computes a discretised approximation to
  the spatial covariance function, using the Fast Fourier Transform.
  The return value is
  another pixel image (object of class \code{"im"}) whose greyscale values
  are values of the spatial covariance function.

  If the argument \code{Y} is present, then \code{imcov(X,Y)}
  computes the set \emph{cross-covariance} function \eqn{C(u)}
  defined as 
  \deqn{
    C(u) = \int X(u)Y(u+v)\, {\rm d}u.
  }{
    C(u) = integral of X(u) * Y(u+v) du.
  }
}
\seealso{
  \code{\link{setcov}},
  \code{\link{owin}},
  \code{\link{as.owin}},
  \code{\link{erosion}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  X <- as.im(w)
  v <- imcov(X)
  plot(v)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
