\name{im}
\alias{im}
\title{Create a Pixel Image Object}
\description{
  Creates an object of
  class \code{"im"} representing a two-dimensional pixel image.
}
\usage{
  im(mat, xcol=seq(ncol(mat)), yrow=seq(nrow(mat)), lev=levels(mat))
}
\arguments{
  \item{mat}{
    matrix or vector containing the pixel values of the image.
  }
  \item{xcol}{
    vector of \eqn{x} coordinates for the pixels
  }
  \item{yrow}{
    vector of \eqn{y} coordinates for the pixels
  }
  \item{lev}{
    possible factor levels, if \code{mat} should be interpreted
    as a factor.
  }
}
\details{
  This function creates an object of class \code{"im"} representing
  a two-dimensional pixel image. See \code{\link{im.object}}
  for details of this class.

  If \code{mat} is not a matrix, it will first be converted into
  a matrix with \code{nrow(mat) = length(yrow)} and
  \code{ncol(mat) = length(xcol)}.
  
  Note that the \bold{row} index of the matrix \code{mat} corresponds
  to increasing \bold{y} coordinate, while the column index of \code{mat}
  corresponds to increasing \bold{x} coordinate. Thus
  \code{xrow} has one entry for each row of \code{mat}
  and \code{ycol} has one entry for each column of \code{mat}. A correct
  display of the image would be obtained by transposing, e.g.
  \code{image.default(xcol, yrow, t(mat))}, if you wanted to do it by hand.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{Kmeasure}},
  \code{\link{setcov}}
}
\section{Warnings}{
  The internal representation of images is likely to change in the
  next release of this package.
}
\examples{
   whitenoise <- im(matrix(rnorm(10000), 100, 100))
   image(whitenoise)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
