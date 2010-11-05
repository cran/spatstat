\name{smooth.ppp}
\alias{smooth.ppp}
\alias{markmean}
\alias{markvar}
\title{Spatial smoothing of observations at irregular points}
\description{
  Performs spatial smoothing of numeric values observed
  at a set of irregular locations. 
}
\usage{
smooth.ppp(X, ..., weights = rep(1, X$n), at="pixels")
markmean(X, ...)
markvar(X, ...)
}
\arguments{
  \item{X}{A marked point pattern (object of class \code{"ppp"}).}
  \item{\dots}{Arguments passed to \code{\link{density.ppp}}
    to control the kernel smoothing and
    the pixel resolution of the result.}
  \item{weights}{Optional weights attached to the observations.}
  \item{at}{
    String specifying whether to compute the intensity values
    at a grid of pixel locations (\code{at="pixels"}) or
    only at the points of \code{X} (\code{at="points"}).
  }
}
\details{
  The function \code{smooth.ppp} performs spatial smoothing of numeric values
  observed at a set of irregular locations. The functions
  \code{markmean} and \code{markvar} are wrappers for \code{smooth.ppp}
  which compute the spatially-varying mean and variance of the marks of
  a point pattern.

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
  of class \code{"ppp"}, see \code{\link{ppp.object}}).
  The points of the pattern are taken to be the
  observation locations \eqn{x_i}{x[i]}, and the marks of the pattern
  are taken to be the numeric values \eqn{v_i}{v[i]} observed at these
  locations.

  The marks are allowed to be a data frame (in \code{smooth.ppp}
  and \code{markmean}). Then the smoothing procedure is applied to each
  column of marks. 
  
  The numerator and denominator are computed by \code{\link{density.ppp}}.
  The arguments \code{...} control the smoothing kernel parameters
  and determine whether edge correction is applied.
  See \code{\link{density.ppp}}.

  The optional argument \code{weights} allows numerical weights to
  be applied to the data. If a weight \eqn{w_i}{w[i]}
  is associated with location \eqn{x_i}{x[i]}, then the smoothed
  function is 
  (ignoring edge corrections)
  \deqn{
    g(u) = \frac{\sum_i k(u-x_i) v_i w_i}{\sum_i k(u-x_i) w_i}
  }{
    g(u) = (sum of k(u-x[i]) v[i] w[i])/(sum of k(u-x[i]) w[i])
  }
}
\value{
  \emph{If \code{X} has a single column of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is
    a pixel image (object of class \code{"im"}). 
    Pixel values are values of the interpolated function.
    \item
    If \code{at="points"}, the result is a numeric vector
    of length equal to the number of points in \code{X}.
    Entries are values of the interpolated function at the points of \code{X}.
  }
  \emph{If \code{X} has a data frame of marks:}
  \itemize{
    \item 
    If \code{at="pixels"} (the default), the result is a named list of 
    pixel images (object of class \code{"im"}). There is one
    image for each column of marks. This list also belongs to
    the class \code{listof}, for which there is a plot method.
    \item
    If \code{at="points"}, the result is a data frame
    with one row for each point of \code{X},
    and one column for each column of marks. 
    Entries are values of the interpolated function at the points of \code{X}.
  }
  The return value has attributes
  \code{"sigma"} and \code{"varcov"} which report the smoothing
  bandwidth that was used.
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
   # mark variance
   plot(markvar(longleaf))
   # data frame of marks: trees marked by diameter and height
   data(finpines)
   plot(smooth.ppp(finpines))
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
