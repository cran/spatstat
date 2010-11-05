\name{Kmeasure}
\alias{Kmeasure}
\title{Reduced Second Moment Measure}
\description{
Estimates the reduced second moment measure \eqn{\kappa}{Kappa} 
from a point pattern in a window of arbitrary shape.
}
\usage{
  Kmeasure(X, sigma, edge=TRUE, ..., varcov=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{\kappa}{Kappa} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
    }
    \item{sigma}{
      Standard deviation \eqn{\sigma}{sigma} of the Gaussian
      smoothing kernel. Incompatible with \code{varcov}.
    }
    \item{edge}{
      logical value indicating whether an edge correction
      should be applied.
    }
    \item{\dots}{Ignored.}
    \item{varcov}{
      Variance-covariance matrix of the Gaussian smoothing kernel.
      Incompatible with \code{sigma}.
    }
}
\value{
  A real-valued pixel image (an object of class \code{"im"},
  see \code{\link{im.object}}) whose pixel values are estimates
  of the value of the reduced second moment measure for each pixel
  (i.e. estimates of the integral of the second moment density
  over each pixel).
}
\details{
  The reduced second moment measure \eqn{\kappa}{Kappa}
  of a stationary point process \eqn{X} is defined so that,
  for a `typical' point \eqn{x} of the process, 
  the expected number of other points \eqn{y} of the process
  such that the vector \eqn{y - x} lies in a region \eqn{A},
  equals \eqn{\lambda \kappa(A)}{lambda Kappa(A)}.
  Here \eqn{\lambda}{lambda}
  is the intensity of the process,
  i.e. the expected number of points of \eqn{X} per unit area.

  The more familiar
  K-function \eqn{K(t)} is just the value of the reduced second moment measure
  for each disc centred at the origin; that is,
  \eqn{K(t) = \kappa(b(0,t))}{K(t) = Kappa(b(0,t))}.

  An estimate of \eqn{\kappa}{Kappa} derived from a spatial point
  pattern dataset can be useful in exploratory data analysis.
  Its advantage over the K-function is that it is also sensitive
  to anisotropy and directional effects. 

  This function computes an estimate of \eqn{\kappa}{Kappa}
  from a point pattern dataset \code{X},
  which is assumed to be a realisation of a stationary point process,
  observed inside a known, bounded window. Marks are ignored.

  The algorithm approximates the point pattern and its window by binary pixel
  images, introduces a Gaussian smoothing kernel
  and uses the Fast Fourier Transform \code{\link{fft}}
  to form a density estimate of \eqn{\kappa}{Kappa}. The calculation
  corresponds to the edge correction known as the ``translation
  correction''.

  The Gaussian smoothing kernel may be specified by either of the
  arguments \code{sigma} or \code{varcov}. If \code{sigma} is a single
  number, this specifies an isotropic Gaussian kernel
  with standard deviation \code{sigma} on each coordinate axis.
  If \code{sigma} is a vector of two numbers, this specifies a Gaussian
  kernel with standard deviation \code{sigma[1]} on the \eqn{x} axis,
  standard deviation \code{sigma[2]} on the \eqn{y} axis, and zero
  correlation between the \eqn{x} and \eqn{y} axes. If \code{varcov} is
  given, this specifies the variance-covariance matrix of the
  Gaussian kernel. There do not seem to be any well-established rules
  for selecting the smoothing kernel in this context.
  
  The density estimate of \eqn{\kappa}{Kappa}
  is returned in the form of a real-valued pixel image.
  Pixel values are estimates of the
  integral of the second moment density over the pixel.
  (The uniform Poisson process would have values identically equal to
  \eqn{a} where \eqn{a} is the area of a pixel.)
  Sums of pixel values over a desired region \eqn{A} are estimates of the
  value of \eqn{\kappa(A)}{Kappa(A)}. The image \code{x} and \code{y}
  coordinates are on the same scale as vector displacements in the
  original point pattern window. The point \code{x=0, y=0} corresponds
  to the `typical point'.
  A peak in the image near \code{(0,0)} suggests clustering;
  a dip in the image near \code{(0,0)} suggests inhibition;
  peaks or dips at other positions suggest possible periodicity.
}
\references{
  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\seealso{
  \code{\link{Kest}},
  \code{\link{fryplot}},
  \code{\link{spatstat.options}},
  \code{\link{im.object}}
}
\examples{
 data(cells)
 image(Kmeasure(cells, 0.05))
 # shows pronounced dip around origin consistent with strong inhibition
 data(redwood)
 image(Kmeasure(redwood, 0.03), col=grey(seq(1,0,length=32)))
 # shows peaks at several places, reflecting clustering and ?periodicity
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
 
 
