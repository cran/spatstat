\name{rescale}
\alias{rescale}
\title{Convert dataset to another unit of length}
\description{
  Converts between different units of length
  in a spatial dataset, such as a point pattern or a window. 
}
\usage{
  rescale(X, s)
}
\arguments{
  \item{X}{Any suitable dataset representing a two-dimensional
    object, such as a point pattern (object of class \code{"ppp"}),
    or a window (object of class \code{"owin"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
}
\value{
  Another object of the same type, representing the same
  data, but expressed in the new units.
}
\details{
  This is generic. Methods are provided for
  point patterns (\code{\link{rescale.ppp}})
  and windows (\code{\link{rescale.owin}}).

  The spatial coordinates in the dataset \code{X} will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.
}
\section{Note}{
  The result of this operation is equivalent to the original dataset.
  If you want to actually change the coordinates by
  a linear transformation, producing a dataset that is not equivalent
  to the original one, use \code{\link{affine}}.
}
\seealso{
  \code{\link{units}},
  \code{\link{rescale.ppp}},
  \code{\link{rescale.owin}},
  \code{\link{affine}},
  \code{\link{rotate}},
  \code{\link{shift}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{math}
