\name{hist.im}
\alias{hist.im}
\title{Histogram of Pixel Values in an Image}
\description{
  Computes and displays a histogram of the pixel values in a pixel image.
  The \code{hist} method for class \code{"im"}.
}
\usage{
  \method{hist}{im}(x, \dots)
}
\arguments{
  \item{x}{A pixel image (object of class \code{"im"}).}
  \item{\dots}{Arguments passed to \code{\link{hist.default}}.}
  }
}
\details{
  This function computes and (by default) displays a histogram
  of the pixel values in the image \code{x}.

  An object of class \code{"im"}
  describes a pixel image. See \code{\link{im.object}})
  for details of this class.

  The function \code{hist.im} is a method for the generic
  function \code{\link{hist}} for the class \code{"im"}. 
  
  Any arguments in \code{...} are passed to \code{\link{hist.default}}.
  For example, such arguments control the axes, and may be used to
  suppress the plotting.
}
\value{
  A single number.
}
\seealso{
  \code{\link{hist}},
  \code{\link{hist.default}},
  \code{\link{im.object}},
  \code{\link{summary.im}}.
}
\examples{
  X <- as.im(function(x,y) {x^2}, unit.square())
  hist(X)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
