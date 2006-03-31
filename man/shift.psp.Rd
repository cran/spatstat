\name{shift.psp}
\alias{shift.psp}
\title{Apply Vector Translation To Line Segment Pattern}
\description{
  Applies a vector shift to a line segment pattern. 
}
\synopsis{
 shift.psp(X, \dots)
}
\usage{
 shift.psp(X, vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Line Segment pattern (object of class \code{"psp"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another line segment pattern (of class \code{"psp"}) representing the
  result of applying the vector shift.
}
\details{
  The line segment pattern, and its window, are
  translated by the vector \code{vec}.
  
  This is a method for the generic function \code{\link{shift}}.
}
\seealso{
  \code{\link{shift}},
  \code{\link{shift.owin}},
  \code{\link{shift.ppp}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  plot(X, col="red")
  Y <- shift(X, c(0.05,0.05))
  plot(Y, add=TRUE, col="blue")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
