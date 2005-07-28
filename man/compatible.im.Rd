\name{compatible.im}
\alias{compatible.im}
\title{Test Whether Two Pixel Images Are Compatible}
\description{
  Tests whether two pixel image objects have compatible dimensions.
}
\usage{
compatible.im(A, B)
}
\arguments{
  \item{A,B}{Pixel images (objects of class \code{"im"}).}
}
\details{
  This function tests whether the pixel images \code{A} and \code{B}
  have compatible pixel dimensions. They are compatible if they have
  the same number of rows and columns, the same physical pixel
  dimensions, and occupy the same rectangle in the plane.
}
\value{
  Logical value: \code{TRUE} if the images are compatible, and \code{FALSE}
  if they are not.
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
