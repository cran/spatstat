\name{as.im}
\alias{as.im}
\title{Convert to Pixel Image}
\description{
  Converts various kinds of data to a pixel image
}
\usage{
  as.im(X, W=as.mask(as.owin(X), dimyx=dimyx), \dots, dimyx=NULL)
}
\arguments{
  \item{X}{Data to be converted to a pixel image.}
  \item{W}{Window object which determines the spatial domain of the image}.
  \item{\dots}{Additional arguments passed to \code{X} when \code{X}
    is a function.}
  \item{dimyx}{
    Optional. Dimensions for the pixel array.
    An integer, or vector of 2 integers \code{c(ny, nx)}.
  }
}
\details{
  This function converts the data \code{X} into a pixel image
  object of class \code{"im"} (see \code{\link{im.object}}).

  Currently \code{X} may be any of the following:
  \itemize{
    \item
    a pixel image object, of class \code{"im"}.
    \item
    a window object, of class \code{"owin"} (see
    \code{\link{owin.object}}). The result is an image
    with all pixel entries equal to \code{1} inside the window \code{X},
    and \code{NA} outside.
    \item
    a single number (or a single logical, complex, factor or character
    value). The result is an image
    with all pixel entries equal to this constant value
    inside the window \code{W} (and \code{NA} outside).
    Argument \code{W} is required.
    \item
    a function of the form \code{function(x, y, ...)}
    which is to be evaluated to yield the image pixel values.
    In this case, the additional argument \code{W} must be present.
    This window will be converted to
    a binary image mask. Then the function \code{X} will be evaluated
    in the form \code{X(x, y, ...)} where \code{x} and \code{y} are
    \bold{vectors} containing the \eqn{x} and \eqn{y} coordinates
    of all the pixels in the image mask, and \code{...} are any extra
    arguments given. This function must return a
    vector or factor of the same length as the input vectors,
    giving the pixel values. 
    \item
    a list with entries \code{x, y, z} in the format expected by
    the standard \code{R} functions
    \code{\link{image.default}} and \code{\link{contour.default}}.
    That is, \code{z} is a matrix of pixel values, \code{x} and \code{y}
    are vectors of \eqn{x} and \eqn{y} coordinates respectively,
    and \code{z[i,j]} is the pixel value for the location
    \code{(x[i],y[j])}. 
  }

  The spatial domain (enclosing rectangle) of the pixel image
  is determined by the argument \code{W}. If \code{W} is absent, 
  the spatial domain is determined by \code{X}.
  When \code{X} is a function or a single numerical value,
  \code{W} is required.
  
  The pixel array dimensions of the final resulting image are determined
  by (in priority order)
  \itemize{
    \item the argument \code{dimyx} if present;
    \item the pixel dimensions of the window \code{W}, if it is
    present and if it is a binary mask;
    \item the pixel dimensions of \code{X} if it is an image,
    a binary mask, or a \code{list(x,y,z)};
    \item the default pixel dimensions, 
    controlled by \code{\link{spatstat.options}}.
  }

  Note that if \code{dimyx} is given, this will override
  the pixel dimensions of \code{X} if it has them.
  Thus, \code{as.im} can be used to change an image's pixel dimensions.
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
  # usual convention in S
  stuff <- list(x=1:10, y=1:10, z=matrix(1:100, nrow=10))
  Z <- as.im(stuff)
  # convert to finer grid
  Z <- as.im(Z, dimyx=256)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{manip}
