\name{print.im}
\alias{print.im}
\title{Print Brief Details of an Image}
\description{
  Prints a very brief description of a pixel image object.
}
\synopsis{
  print.im(x, \dots)
}
\usage{
 print.im(x)
 x
}
\arguments{
  \item{x}{Pixel image (object of class \code{"im"}).}
}
\details{
  A very brief description of the pixel image \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{im.object}},
  \code{\link{summary.im}}
}
\examples{
  data(letterR)
  U <- as.im(letterR)
  U
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
