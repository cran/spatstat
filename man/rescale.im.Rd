\name{rescale.im}
\alias{rescale.im}
\title{Convert Pixel Image to Another Unit of Length}
\description{
  Converts a pixel image to 
  another unit of length.
}
\usage{
 \method{rescale}{im}(X, s)
}
\arguments{
  \item{X}{Pixel image (object of class \code{"im"}).}
  \item{s}{Conversion factor: the new units are \code{s} times the old units.}
}
\value{
  Another pixel image (of class \code{"im"}),
  containing the same pixel values, but with pixel coordinates
  expressed in the new units.
}
\details{
  This is a method for the generic function \code{\link{rescale}}.

  The spatial coordinates of the pixels in \code{X}
  will be re-expressed
  in terms of a new unit of length that is \code{s} times the current
  unit of length given in \code{X}.
  (Thus, the coordinate values are \emph{divided} by \code{s},
  while the unit value is multiplied by \code{s}).

  The result is a pixel image representing the \emph{same} data
  but re-expressed in a different unit.

  Pixel values are unchanged. 
}
\seealso{
  \code{\link{im}},
  \code{\link{rescale}},
  \code{\link{unitname}},
  \code{\link{eval.im}}
}
\examples{
# Bramble Canes data: 1 unit = 9 metres
  data(bramblecanes)
# distance transform
  Z <- distmap(bramblecanes)
# convert to metres
# first alter the pixel values
  Zm <- eval.im(9 * Z)
# now rescale the pixel coordinates
  Zm <- rescale(Zm, 1/9)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
