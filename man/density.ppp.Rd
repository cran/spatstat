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
    Point pattern (object of class \code{"ppp"}).
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
  Pixel values are estimated intensity values,
  expressed in \dQuote{points per unit area}.
}
\details{
  This is a method for the generic function \code{density}.

  It computes a fixed-bandwidth kernel estimate 
  (Diggle, 1985) of the intensity function of the point process
  that generated the point pattern \code{x}.
  
  By default it computes the convolution of the
  isotropic Gaussian kernel of standard deviation \code{sigma}
  with point masses at each of the data points in \code{x}.
  Each point has unit weight, unless the argument \code{weights} is
  given (it should be a numeric vector; weights can be negative or zero).

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

  To perform spatial interpolation of values that were observed
  at the points of a point pattern, use \code{\link{smooth.ppp}}.

  For adaptive nonparametric estimation, see \code{\link{adaptive.density}}.
}
\seealso{
  \code{\link{smooth.ppp}},
  \code{\link{adaptive.density}},
  \code{\link{ppp.object}},
  \code{\link{im.object}}
}
\note{
  This function is often misunderstood.

  The result of \code{density.ppp} is not a spatial smoothing 
  of the marks or weights attached to the point pattern.
  To perform spatial interpolation of values that were observed
  at the points of a point pattern, use \code{\link{smooth.ppp}}.

  The result of \code{density.ppp} is not a probability density.
  It is an estimate of the \emph{intensity function} of the
  point process that generated the point pattern data.
  Intensity is the expected number of random points
  per unit area.
  The units of intensity are \dQuote{points per unit area}.
  Intensity is usually a function of spatial location,
  and it is this function which is estimated by \code{density.ppp}.
  The integral of the intensity function over a spatial region gives the
  expected number of points falling in this region.

  Inspecting an estimate of the intensity function is usually the
  first step in exploring a spatial point pattern dataset.
  For more explanation, see the workshop notes (Baddeley, 2008)
  or Diggle (2003).
}
\examples{
  data(cells)
  Z <- density.ppp(cells, 0.05)
  plot(Z)
}
\references{
  Baddeley, A. (2008) Analysing spatial point patterns in R.
  Workshop notes. CSIRO online technical publication.
  URL: \code{www.csiro.au/resources/pf16h.html}

  Diggle, P.J. (1985)
  A kernel method for smoothing point process data.
  \emph{Applied Statistics} (Journal of the Royal Statistical Society,
  Series C) \bold{34} (1985) 138--147.

  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
