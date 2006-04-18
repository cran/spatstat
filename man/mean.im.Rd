\name{mean.im}
\alias{mean.im}
\title{Mean Pixel Value in an Image}
\description{
  Calculates the mean of the pixel values in a pixel image.
  The \code{mean} method for class \code{"im"}.
}
\usage{
  \method{mean}{im}(x, \dots)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{Arguments passed to \code{\link{mean.default}}.}
  }
}
\details{
  This function calculates the mean value of the pixels in the image
  \code{x}.

  An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class.

  The function \code{mean.im} is a method for the generic
  function \code{\link{mean}} for the class \code{"im"}. 
  
  If the image \code{x} is logical-valued, the mean value of \code{x} is
  the fraction of pixels that have the value \code{TRUE}.
  If the image \code{x} is factor-valued, then the mean of \code{x}
  is the mean of the integer codes of the pixel values.

  Any arguments in \code{...} are passed to \code{\link{mean.default}}.
  In particular, using the argument \code{trim} will compute the
  trimmed mean, as explained in the help for \code{\link{mean.default}}.

  Other information about an image can be obtained using
  \code{\link{summary.im}}.
}
\value{
  A single number.
}
\seealso{
  \code{\link{mean}},
  \code{\link{mean.default}},
  \code{\link{im.object}},
  \code{\link{summary.im}}.
}
\examples{
  X <- as.im(function(x,y) {x^2}, unit.square())
  mean(X)
  mean(X, trim=0.05)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
