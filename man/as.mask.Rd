\name{as.mask}
\alias{as.mask}
\title{Pixel Image Approximation of a Window}
\description{
  Obtain a discrete (pixel image) approximation of a given window
}
\usage{
 as.mask(w, eps=NULL, dimyx=NULL)
}
\arguments{
  \item{w}{A window (object of class \code{"owin"}).}
  \item{eps}{(optional) width and height of pixels.}
  \item{dimyx}{(optional) pixel array dimensions}
}
\value{
  A window (object of class \code{"owin"})
  of type \code{"mask"} representing a binary pixel image.
}
\details{
  The shape of the window \code{w} will be approximated
  by a binary pixel image. See \code{\link{owin.object}} for
  explanation of binary image masks.

  The pixel raster will be an \eqn{m \times n}{m x n} rectangular grid
  where \eqn{m, n} are given by \code{dimyx[2]}, \code{dimyx[1]}
  respectively (if the argument \code{dimyx} is given)
  or by the constraint that
  pixel widths and heights should be approximately \code{eps}.
  If neither \code{eps} nor \code{dimyx} is given,
  the pixel raster dimensions are obtained from
  \code{spatstat.options("npixel")}.

  There is no inverse of this function!
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{as.rectangle}},
  \code{\link{spatstat.options}}
}
\examples{
  \dontrun{library(spatstat)}
  w <- owin(c(0,10),c(0,10), poly=list(x=c(1,2,3,2,1), y=c(2,3,4,6,7)))
  \dontrun{plot(w)}
  m <- as.mask(w)
  \dontrun{plot(m)}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
