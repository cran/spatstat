\name{shift.im}
\alias{shift.im}
\title{Apply Vector Translation To Pixel Image}
\description{
  Applies a vector shift to a pixel image
}
\usage{
 \method{shift}{im}(X, vec=c(0,0), \dots)
}
\arguments{
  \item{X}{Pixel image (object of class \code{"im"}).}
  \item{vec}{Vector of length 2 representing a translation.}
  \item{\dots}{Ignored}
}
\value{
  Another pixel image (of class \code{"im"}) representing the
  result of applying the vector shift.
}
\details{
  The spatial location of each pixel in the image
  is translated by the vector \code{vec}.
  This is a method for the generic function \code{\link{shift}}.
}
\seealso{
  \code{\link{shift}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 Y <- shift(X, c(10,10))
 plot(Y)
 # no discernible difference except coordinates are different
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
