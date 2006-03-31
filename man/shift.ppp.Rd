\name{shift.ppp}
\alias{shift.ppp}
\title{Apply Vector Translation To Point Pattern}
\description{
  Applies a vector shift to a point pattern. 
}
\usage{
 shift.ppp(X, vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another point pattern (of class \code{"ppp"}) representing the
  result of applying the vector shift.
}
\details{
  The point pattern, and its window, are
  translated by the vector \code{vec}.
  
  This is a method for the generic function \code{\link{shift}}.
}
\seealso{
  \code{\link{shift}},
  \code{\link{shift.owin}},
  \code{\link{rotate}},
  \code{\link{affine}}
}
\examples{
  data(cells)
  X <- shift(cells, c(2,3))
  \dontrun{
  plot(X)
  # no discernible difference except coordinates are different
  }
  Y <- superimpose(cells, shift(cells, c(0.1,0.1)))
  \dontrun{plot(Y)}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
