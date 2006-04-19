\name{affine.owin}
\alias{affine.owin}
\title{Apply Affine Transformation To Window}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a window. 
}
\usage{
  \method{affine}{owin}(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Window (object of class \code{"owin"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another window (of class \code{"owin"}) representing the
  result of applying the affine transformation.
}
\details{
  The window is subjected first to the linear transformation represented by
  \code{mat} (multiplying on the left by \code{mat}),
  and then the result is translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix.

  This is a method for the generic function \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
  \code{\link{affine.ppp}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  # shear transformation
  X <- affine(owin(), matrix(c(1,0,0.6,1),ncol=2))
  \dontrun{
  plot(X)
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
