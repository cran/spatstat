\name{density.ppp}
\alias{density.ppp}
\title{Kernel Smoothed Intensity of Point Pattern}
\description{
  Compute a kernel smoothed intensity function from a point pattern.
}
\usage{
  \method{density}{ppp}(x, sigma, \dots,
        weights, edge=TRUE, varcov=NULL,
        at="pixels", leaveoneout=TRUE)
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
  \item{at}{
    String specifying whether to compute the intensity values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{x} (\code{at="points"}).
  }
  \item{leaveoneout}{
    Logical value indicating whether to compute a leave-one-out
    estimator. Applicable only when \code{at="points"}.
  }
}
\value{
  By default, the result is
  a pixel image (object of class \code{"im"}). 
  Pixel values are estimated intensity values,
  expressed in \dQuote{points per unit area}.

  If \code{at="points"}, the result is a numeric vector
  of length equal to the number of points in \code{x}.
  Values are estimated intensity values at the points of \code{x}.

  In either case, the return value has attributes
  \code{"sigma"} and \code{"varcov"} which report the smoothing
  bandwidth that was used.
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

  Thus the intensity value at a point \eqn{u} is
  \deqn{
    \hat\lambda(u) = e(u) \sum_i k(x_i - u) w_i
  }{
    lambda(u) = e(u) sum[i] k(x[i] - u) w[i]
  }
  where \eqn{k} is the Gaussian smoothing kernel,
  \eqn{e(u)} is an edge correction factor, 
  and \eqn{x_i}{w[i]} are the weights.

  By default the intensity values are
  computed at every location \eqn{u} in a fine grid,
  and are returned as a pixel image.
  Computation is performed using the Fast Fourier Transform.
  Accuracy depends on the pixel resolution, controlled by the arguments
  \code{\dots} passed to \code{\link{as.mask}}.

  If \code{at="points"}, the intensity values are computed 
  to high accuracy at the points of \code{x} only. Computation is
  performed by directly evaluating and summing the Gaussian kernel
  contributions without discretising the data. The result is a numeric
  vector giving the density values.
  The intensity value at a point \eqn{x_i}{x[i]} is
  \deqn{
    \hat\lambda(x_i) = e(x_i) \sum_j k(x_j - x_i) w_i
  }{
    lambda(x[i]) = e(x[i]) sum[j] k(x[j] - x[i]) w[i]
  }
  If \code{leaveoneout=TRUE} (the default), then the sum in the equation
  is taken over all \eqn{j} not equal to \eqn{i},
  so that the intensity value at a
  data point is the sum of kernel contributions from
  all \emph{other} data points.
  If \code{leaveoneout=FALSE} then the sum is taken over all \eqn{j},
  so that the intensity value at a data point includes a contribution
  from the same point.

  To perform spatial interpolation of values that were observed
  at the points of a point pattern, use \code{\link{smooth.ppp}}.

  For adaptive nonparametric estimation, see
  \code{\link{adaptive.density}}.
  For data sharpening, see \code{\link{sharpen.ppp}}.
}
\seealso{
  \code{\link{smooth.ppp}},
  \code{\link{sharpen.ppp}},
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
  Z <- density(cells, 0.05)
  plot(Z)
  density(cells, 0.05, at="points")
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
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
