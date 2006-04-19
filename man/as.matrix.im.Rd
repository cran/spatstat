\name{as.matrix.im}
\alias{as.matrix.im}
\title{Convert Pixel Image to Matrix}
\description{
  Convert a pixel image to a matrix
}
\usage{
  \method{as.matrix}{im}(x)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
}
\details{
  This function takes the pixel image \code{x}
  and returns a matrix containing the pixel values.

  It is handy when you want to extract a summary of the pixel values.
  See the Examples.
}
\value{
  A matrix.
}
\examples{
   # artificial image
   Z <- setcov(square(1))

   M <- as.matrix(Z)

   median(M)
   
   \dontrun{
   # plot the cumulative distribution function of pixel values
   plot(ecdf(as.matrix(Z)))
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
