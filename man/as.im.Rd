\name{as.im}
\alias{as.im}
\title{Convert to Pixel Image}
\description{
  Converts various kinds of data to a pixel image
}
\usage{
  as.im(X, W, \dots)
}
\arguments{
  \item{X}{Data to be converted to a pixel image.}
  \item{W}{Window object required when \code{X} is a function.}
  \item{\dots}{Additional arguments passed to \code{X} when \code{X}
    is a function.}
}
\details{
  This function converts the data \code{x} into a pixel image
  object of class \code{"im"} (see \code{\link{im.object}}).

  Currently \code{X} may be any of the following:
  \itemize{
    \item
    a pixel image object, of class \code{"im"}. This object is just returned.
    \item
    a window object, of class \code{"owin"} (see \code{\link{owin.object}}).
    This window is first converted to a binary image mask
    by \code{\link{as.mask}}. Then the pixel entries 
    \code{TRUE} and \code{FALSE} are converted to
    \code{1} and \code{NA} respectively.
    The result is returned as an object of class \code{"im"}.
    \item
    a function of the form \code{function(x, y, ...)}
    which is to be evaluated to yield the image pixel values.
    In this case, the additional argument \code{W} must be
    a window object. This window will be converted to
    a binary image mask. Then the function \code{X} will be evaluated
    in the form \code{X(x, y, ...)} where \code{x} and \code{y} are
    \bold{vectors} containing the \eqn{x} and \eqn{y} coordinates
    of all the pixels in the image mask. This must return a numeric vector
    of the same length as the input vectors, giving the pixel values.
    \item
    a single numerical value. The result is an image
    with pixel entries equal to this numerical value
    inside the image \code{win}.
  }
}
\value{
  An image object of class \code{"im"}.
}
\examples{
  data(demopat)
  # window object
  W <- demopat$window
  plot(W)
  Z <- as.im(W)
  image(Z)
  # function
  Z <- as.im(function(x,y) {x^2 + y^2}, unit.square())
  image(Z)
  # function with extra arguments
  f <- function(x, y, x0, y0) {
      sqrt((x - x0)^2 + (y-y0)^2)
  }
  Z <- as.im(f, unit.square(), x0=0.5, y0=0.5)
  image(Z)
  # Revisit the Sixties
  data(letterR)
  Z <- as.im(f, letterR, x0=2.5, y0=2)
  image(Z)  
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
