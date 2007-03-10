\name{subsetget.im}
\alias{[<-.im}
\title{Reset Values in Subset of Image}
\description{
  Reset the values in a subset of a pixel image.
}
\usage{
  \method{[}{im}(x, i) <- value
}
\arguments{
  \item{x}{
    A two-dimensional pixel image.
    An object of class \code{"im"}.
  }
  \item{i}{
    Object defining the subregion or subset to be extracted.
    Either a spatial window (an object of class \code{"owin"})
    or a point pattern (an object of class \code{"ppp"}).
  }
  \item{value}{
     Vector of replacement values (short vectors will be
     recycled).
  }
} 
\value{
  The image \code{x} with the values replaced.
}
\details{
  This function changes some of the pixel values in a
  pixel image. The image \code{X} must be an object of class
  \code{"im"} representing a pixel image defined inside a
  rectangle in two-dimensional space (see \code{\link{im.object}}).

  The subset to be changed is determined by the argument \code{i}.
  If \code{i} is a spatial window (an object of class \code{"owin"}),
  the values of the image inside this window are changed.
  If \code{i} is a point pattern (an object of class
  \code{"ppp"}), then the values of the pixel image at the points of
  this pattern are changed.
}
\seealso{
  \code{\link{im.object}},
  \code{\link{[.im}},
  \code{\link{ppp.object}},
  \code{\link{owin.object}}
}
\examples{
 # make up an image
 X <- setcov(unit.square())
 plot(X)

 # a rectangular subset
 W <- owin(c(0,0.5),c(0.2,0.8))
 X[W] <- 2
 plot(X)

 # a polygonal subset
 data(letterR)
 R <- affine(letterR, diag(c(1,1)/2), c(-2,-0.7))
 X[R] <- 3
 plot(X)

 # a point pattern
 P <- rpoispp(20)
 X[P] <- 10
 plot(X)
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
