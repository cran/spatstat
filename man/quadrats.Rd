\name{quadrats}
\alias{quadrats}
\title{Divide Region into Quadrats}
\description{
  Divides window into rectangular quadrats and 
  returns the quadrats as a tessellation.
}
\usage{
quadrats(X, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL)
}
\arguments{
  \item{X}{
    A window (object of class \code{"owin"})
    or anything that can be coerced to a window using
    \code{\link{as.owin}}, such as a point pattern.
  }
  \item{nx,ny}{
    Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
    Incompatible with \code{xbreaks} and \code{ybreaks}.
  }
  \item{xbreaks}{
    Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the quadrats. Incompatible with \code{ny}.
  }
}
\details{
  The window \code{X} is divided into
  an \code{nx * ny} grid of rectangular tiles or `quadrats'.
  These tiles are returned as a tessellation (object of class
  \code{"tess"}) which can be plotted and used in other analyses.

  If \code{xbreaks} is given, it should be a numeric vector
  giving the \eqn{x} coordinates of the quadrat boundaries.
  If it is not given, it defaults to a
  sequence of \code{nx+1} values equally spaced
  over the range of \eqn{x} coordinates in the window \code{X$window}.

  Similarly if \code{ybreaks} is given, it should be a numeric
  vector giving the \eqn{y} coordinates of the quadrat boundaries.
  It defaults to a vector of \code{ny+1} values
  equally spaced over the range of \eqn{y} coordinates in the window.
  The lengths of \code{xbreaks} and \code{ybreaks} may be different.
}
\value{
  A tessellation (object of class \code{"tess"}) as described under
  \code{\link{tess}}.
}
\examples{
 W <- square(10)
 Z <- quadrats(W, 4, 5)
 plot(Z)
}
\seealso{
  \code{\link{tess}},
  \code{\link{quadratcount}},
  \code{\link{quadrat.test}},
  \code{\link{quadratresample}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{utilities}
\keyword{datagen}
