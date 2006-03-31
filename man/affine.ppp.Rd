\name{affine.ppp}
\alias{affine.ppp}
\title{Apply Affine Transformation To Point Pattern}
\description{
  Applies any affine transformation of the plane (linear transformation
  plus vector shift) to a point pattern. 
}
\usage{
 affine.ppp(X, mat=diag(c(1,1)), vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{mat}{Matrix representing a linear transformation.}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another point pattern (of class \code{"ppp"}) representing the
  result of applying the affine transformation.
}
\details{
  The point pattern, and its window, are subjected first to the
  linear transformation represented by
  \code{mat} (multiplying on the left by \code{mat}),
  and are then translated by the vector \code{vec}.
  
  The argument \code{mat} must be a nonsingular \eqn{2 \times 2}{2 * 2}
  matrix.

  This is a method for the generic function \code{\link{affine}}.
}
\seealso{
  \code{\link{affine}},
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
