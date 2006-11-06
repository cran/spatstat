\name{density.ppp}
\alias{density.ppp}
\title{Kernel Smoothed Intensity of Point Pattern}
\description{
  Compute a kernel smoothed intensity function from a point pattern.
}
\usage{
  \method{density}{ppp}(x, sigma, \dots, weights, edge=TRUE, varcov=NULL)
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
  \item{\dots}{
    Arguments passed to \code{\link{as.mask}} to determine
    the pixel resolution.
  }
  \item{edge}{
    Logical flag: if \code{TRUE}, apply edge correction.
  }
  \item{varcov}{
    Variance-covariance matrix of anisotropic Gaussian kernel.
    Incompatible with \code{sigma}.
  }
}
\value{
  A pixel image (object of class \code{"im"}).
}
\details{
  This is a method for the generic function \code{density}.
  
  A kernel estimate of the intensity function of the point pattern
  is computed. The result is 
  the convolution of the isotropic Gaussian kernel of
  standard deviation \code{sigma} with point masses at each of the data
  points. The default is to assign
  a unit weight to each point.
  If \code{weights} is present, the point masses have these
  weights (which may be signed real numbers).

  If \code{edge=TRUE}, the intensity estimate is corrected for
  edge effect bias by dividing it by the convolution of the
  Gaussian kernel with the window of observation.

  Instead of the isotropic Gaussian kernel with standard deviation
  \code{sigma}, the smoothing kernel may be chosen to be any Gaussian
  kernel, by giving the variance-covariance matrix \code{varcov}.
  The arguments \code{sigma} and \code{varcov} are incompatible.
  Also \code{sigma} may be a vector of length 2 giving the
  standard deviations of two independent Gaussian coordinates,
  thus equivalent to \code{varcov = diag(sigma^2)}.
  
  Computation is performed using the Fast Fourier Transform.
  Accuracy depends on the pixel resolution, controlled by the arguments
  \code{\dots} passed to \code{\link{as.mask}}.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\examples{
  data(cells)
  Z <- density.ppp(cells, 0.05)
  plot(Z)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
