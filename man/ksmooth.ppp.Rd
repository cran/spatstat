\name{ksmooth.ppp}
\alias{ksmooth.ppp}
\title{Kernel Smoothed Intensity of Point Pattern}
\description{
  Compute a kernel smoothed intensity function from a point pattern.
}
\synopsis{
  ksmooth.ppp(x, sigma, ..., edge=TRUE)
}
\usage{
  ksmooth.ppp(x, sigma, weights, edge=TRUE)
}
\arguments{
  \item{x}{
    Point pattern (object of class \code{"ppp"}) to be smoothed.
  }
  \item{sigma}{
    Standard deviation of isotropic Gaussian smoothing kernel.
  }
  \item{weights}{
    Optional vector of weights to be attached to the points.
    May include negative values. 
  }
  \item{edge}{
    Logical flag: if \code{TRUE}, apply edge correction.
  }
}
\value{
  A pixel image (object of class \code{"im"}).
}
\details{
  A kernel estimate of the intensity of the point pattern
  is computed. The result is 
  the convolution of the isotropic Gaussian kernel with
  standard deviation \code{sigma} and point masses at each of the data
  points.   The default is to assign
  a unit weight to each point.
  If \code{weights} is present, the point masses have these
  weights (which may be signed real numbers).

  Computation is performed using the Fast Fourier Transform.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\examples{
  data(cells)
  Z <- ksmooth.ppp(cells, 0.05)
  plot(Z)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}}
\keyword{spatial}
