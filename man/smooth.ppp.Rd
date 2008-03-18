\name{smooth.ppp}
\alias{smooth.ppp}
\title{Spatial smoothing of observations at irregular points}
\description{
  Performs spatial smoothing of numeric values observed
  at a set of irregular locations. 
}
\usage{
smooth.ppp(X, ..., weights = rep(1, X$n))
}
\arguments{
  \item{X}{A marked point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Arguments passed to \code{\link{density.ppp}}
    to control the kernel smoothing.}
  \item{weights}{Optional weights attached to the observations.}
}
\details{
  This function performs spatial smoothing of numeric values
  observed at a set of irregular locations.

  Smoothing is performed by Gaussian kernel weighting. If the
  observed values are \eqn{v_1,\ldots,v_n}{v[1],...,v[n]}
  at locations \eqn{x_1,\ldots,x_n}{x[1],...,x[n]} respectively,
  then the smoothed value at a location \eqn{u} is
  (ignoring edge corrections)
  \deqn{
    g(u) = \frac{\sum_i k(u-x_i) v_i}{\sum_i k(u-x_i)}
  }{
    g(u) = (sum of k(u-x[i]) v[i])/(sum of k(u-x[i]))
  }
  where \eqn{k} is a Gaussian kernel.
  
  The argument \code{X} must be a marked point pattern (object
  of class \code{"ppp"}, see \code{\link{ppp.object}})
  in which the points are the observation locations,
  and the marks are the numeric values observed at each point.
  
  The numerator and denominator are computed by \code{\link{density.ppp}}.
  The arguments \code{...} control the smoothing kernel parameters
  and determine whether edge correction is applied.
  See \code{\link{density.ppp}}.

  The optional argument \code{weights} allows numerical weights to
  be applied to the data (the weights appear in both the sums
  in the equation above).
}
\value{
  A pixel image (object of class \code{"im"}).
  Pixel values are values of the interpolated function.
}
\seealso{
  \code{\link{density.ppp}},
  \code{\link{ppp.object}},
  \code{\link{im.object}}.
  To perform interpolation, see the \code{akima} package.
}
\examples{
   # Longleaf data - tree locations, marked by tree diameter
   data(longleaf)
   # Local smoothing of tree diameter
   Z <- smooth.ppp(longleaf)
   # Kernel bandwidth sigma=5
   plot(smooth.ppp(longleaf, 5))
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
