\name{im.object}
\alias{im.object}
\title{Class of Images}
\description{
  A class \code{"im"} to represent a two-dimensional pixel image.
}
\details{
  An object of this class represents
  a two-dimensional pixel image. It specifies
  \itemize{
    \item the dimensions of the rectangular array of pixels
    \item \eqn{x} and \eqn{y} coordinates for the pixels
    \item a numeric value (``grey value'') at each pixel
  }
  If \code{X} is an object of type \code{im},
  it contains the following elements:
  \tabular{ll}{
    \code{v} \tab matrix of values \cr
    \code{dim} \tab dimensions of matrix \code{v} \cr
    \code{xrange} \tab range of \eqn{x} coordinates of image window \cr
    \code{yrange} \tab range of \eqn{y} coordinates of image window \cr
    \code{xstep} \tab width of one pixel \cr
    \code{ystep} \tab height of one pixel \cr
    \code{xcol} \tab vector of \eqn{x} coordinates of centres of pixels \cr
    \code{yrow} \tab vector of \eqn{y} coordinates of centres of pixels 
  }
  Users are strongly advised not to manipulate these entries
  directly.

  Objects of class \code{"im"}
  may be created by the function
  \code{\link{im}}. They are also returned by the functions
  \code{\link{Kmeasure}} and
  \code{\link{setcov}}.
  They may be displayed using the methods
  \code{image.im}, \code{persp.im}
  and \code{contour.im}.
  Values of an image may be looked-up using
  \code{lookup.im}.

  Note that the \bold{row} index of the matrix \code{v} corresponds
  to increasing \bold{y} coordinate, while the column index of \code{m}
  corresponds to increasing \bold{x} coordinate. Thus a correct
  display of the image would be obtained by transposing, e.g.
  \code{image.default(xcol, yrow, t(im))}, if you wanted to do it by hand.
}
\seealso{
  \code{\link{im}},
  \code{\link{Kmeasure}},
  \code{\link{setcov}}
}
\section{Warnings}{
  The internal representation of images is likely to change in the
  next release of this package.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
