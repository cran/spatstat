\name{im}
\alias{im}
\title{Create a Pixel Image Object}
\description{
  Creates an object of
  class \code{"im"} representing a two-dimensional pixel image.
}
\usage{
  im(mat, xcol=seq(ncol(mat)), yrow=seq(nrow(mat)), lev=levels(mat),
   units=NULL)
}
\arguments{
  \item{mat}{
    matrix or vector containing the pixel values of the image.
  }
  \item{xcol}{
    vector of \eqn{x} coordinates for the pixel gid
  }
  \item{yrow}{
    vector of \eqn{y} coordinates for the pixel grid
  }
  \item{lev}{
    possible factor levels, if \code{mat} should be interpreted
    as a factor.
  }
  \item{units}{
    Optional. Name of unit of length. Either a single character string,
    or a vector of two character strings giving the
    singular and plural forms, respectively.
  }
}
\details{
  This function creates an object of class \code{"im"} representing
  a two-dimensional pixel image. See \code{\link{im.object}}
  for details of this class.

  The matrix \code{mat} contains the \sQuote{greyscale} values
  for a rectangular grid of pixels.
  Note carefully that the entry \code{mat[i,j]}
  gives the pixel value at the location \code{(xcol[j],yrow[i])}.
  That is, the \bold{row} index of the matrix \code{mat} corresponds
  to increasing \bold{y} coordinate, while the column index of \code{mat}
  corresponds to increasing \bold{x} coordinate.
  Thus \code{yrow} has one entry for each row of \code{mat}
  and \code{xcol} has one entry for each column of \code{mat}.
  Under the usual convention in \R, a correct
  display of the image would be obtained by transposing the matrix, e.g.
  \code{image.default(xcol, yrow, t(mat))}, if you wanted to do it by hand.

  The entries of \code{mat} may be numeric (real or integer), complex, 
  logical, character, or factor values.
  If \code{mat} is not a matrix, it will be converted into
  a matrix with \code{nrow(mat) = length(yrow)} and
  \code{ncol(mat) = length(xcol)}.

  \R does not allow a matrix to have
  factor-valued entries. So to make a factor-valued image from raw data,
  you must supply \code{mat} as a factor vector and specify the arguments
  \code{xcol} and \code{yrow} to determine the dimensions of the image.
  See the examples. (Alternatively you can use the functions
  \code{\link{cut.im}} or \code{\link{eval.im}} to make factor-valued
  images from other images).

  For a description of the methods available for pixel image objects,
  see \code{\link{im.object}}.

  To convert other kinds of data to a pixel image (for example,
  functions or windows), use \code{\link{as.im}}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{as.im}},
  \code{\link{as.matrix.im}},
  \code{\link{[.im}},
  \code{\link{eval.im}}
}
\section{Warnings}{
  The internal representation of images is likely to change in future
  releases of \pkg{spatstat}. The safe way to extract pixel values
  from an image object is to use \code{\link{as.matrix.im}}
  or \code{\link{[.im}}.
}
\examples{
   vec <- rnorm(1200)
   mat <- matrix(vec, nrow=30, ncol=40)
   whitenoise <- im(mat)
   whitenoise <- im(mat, xcol=seq(0,1,length=40), yrow=seq(0,1,length=30))
   whitenoise <- im(vec, xcol=seq(0,1,length=40), yrow=seq(0,1,length=30))
   plot(whitenoise)

   # FACTOR-VALUED IMAGES:
   cutvec <- cut(mat, 3)
   # although mat was a matrix, cutvec is a vector, with factor values
   cutwhite <- im(cutvec, xcol=seq(0,1,length=40), yrow=seq(0,1,length=30))
   # cutwhite is a factor-valued image
   plot(cutwhite)
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
\keyword{datagen}
 
 
