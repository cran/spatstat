\name{affine}
\alias{affine}
\title{Apply Affine Transformation}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a plane geometrical object,
  such as a point pattern or a window. 
}
\synopsis{
  affine(X, \dots)
}
\usage{
 affine(X, mat=diag(c(1,1)), vec=c(0,0))
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
}
\value{
  Another object of the same type, representing the
  result of applying the affine transformation.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{affine.ppp}})
  and windows (\code{\link{affine.owin}}).

  First the linear transformation represented by
  \code{mat} is applied (multiplying on the left by \code{mat}),
  and then the result is translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix. 
}
\seealso{
  \code{\link{affine.ppp}},
  \code{\link{affine.owin}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\examples{
  data(cells)
  # shear transformation
  X <- affine(cells, matrix(c(1,0,0.6,1),ncol=2))
  \dontrun{
  plot(X)
  # rescale y coordinates by factor 1.3
  plot(affine(cells, diag(c(1,1.3))))
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
