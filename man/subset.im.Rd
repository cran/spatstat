\name{subset.im}
\alias{subset.im}
\alias{[.im}
\title{Extract Subset of Image}
\description{
  Extract a subset or subregion of a pixel image.
}
\synopsis{
  \method{[}{im}(x, i, j, drop=TRUE, \dots)
}
\usage{
  \method{[}{im}(x, w)
  \method{[}{im}(x, w, drop=TRUE)
}
\arguments{
  \item{x}{
    A two-dimensional pixel image.
    An object of class \code{"im"}.
  }
  \item{w}{
    Object defining the subregion or subset to be extracted.
    Either a spatial window (an object of class \code{"owin"})
    or a point pattern (an object of class \code{"ppp"}).
  }
  \item{drop}{
    Logical value. 
    Locations in \code{w} that lie outside the spatial domain of the
    image \code{X} return a pixel value of \code{NA} if
    \code{drop=FALSE}, and are omitted if \code{drop=TRUE}.
  }
} 
\value{
  Either a pixel image or a vector of pixel values. See Details.
}
\details{
  This function extracts a subset of the pixel values in a
  pixel image. The image \code{X} must be an object of class
  \code{"im"} representing a pixel image defined inside a
  rectangle in two-dimensional space (see \code{\link{im.object}}).

  The subset to be extracted is determined by the argument \code{w}.
  If \code{w} is a spatial window (an object of class \code{"owin"}),
  the values of the image inside this window are extracted
  (after first clipping the window to the spatial domain of the image
  if necessary). If \code{w} is a point pattern (an object of class
  \code{"ppp"}), then the values of the pixel image at the points of
  this pattern are extracted.
  
  At locations outside the spatial domain of the image,
  the pixel value is undefined, and is taken to be \code{NA}. 
  The logical argument \code{drop} determines whether such \code{NA}
  values will be returned or omitted.
  It also influences the format of the return value. 

  If \code{w} is a point pattern, 
  \code{X[w, drop=FALSE]} is a numeric vector
  containing the pixel values at each of the points of the pattern.
  Its length is equal to the number of points in the pattern \code{w}.
  It may contain \code{NA}s corresponding to points which lie outside
  the spatial domain of the image \code{X}.
  By contrast, \code{X[w]} or \code{X[w, drop=TRUE]} contains only those
  pixel values which are not \code{NA}. It may be shorter.

  If \code{w} is a spatial window
  then \code{X[w, drop=FALSE]} is another pixel image
  of the same dimensions as \code{X} obtained by 
  setting all pixels outside the window \code{w} to have value
  \code{NA}. When the result is displayed by \code{\link{plot.im}} the effect
  is that the pixel image \code{X} is clipped to the window \code{w}.

  If \code{w} is a spatial window which is \emph{not} a rectangle
  (\code{w$type != "rectangle"}) then \code{X[w, drop=TRUE]}
  is a numeric vector containing the pixel values for all pixels
  that lie inside the window \code{w}.

  If \code{w} is a rectangle (a spatial window
  with \code{w$type = "rectangle"}) then \code{X[w, drop=TRUE]}
  is a pixel image. The spatial domain of this image is the
  intersection of \code{w} with the spatial domain of \code{X}.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{ppp.object}},
  \code{\link{owin.object}},
  \code{\link{plot.im}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 # a rectangular subset
 W <- owin(c(0,0.5),c(0.2,0.8))
 Y <- X[W]
 plot(Y)

 # a polygonal subset
 data(letterR)
 R <- affine(letterR, diag(c(1,1)/2), c(-2,-0.7))
 Y <- X[R, drop=FALSE]
 plot(Y)

 # a point pattern
 P <- rpoispp(20)
 Y <- X[P]
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
