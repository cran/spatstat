\name{Linhom}
\alias{Linhom}
\title{L-function}
\description{
  Calculates an estimate of the inhomogeneous version of Ripley's L-function
  for a spatial point pattern.
}
\usage{
  Linhom(...)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link{Kinhom}}
    to estimate the inhomogeneous K-function.
  }
}
\details{
  This command computes an estimate of the inhomogeneous version of
  the L-function for a spatial point pattern

  The original L-function is a transformation of Ripley's K-function,
  \deqn{L(r) = \sqrt{\frac{K(r)}{\pi}}}{L(r) = sqrt(K(r)/pi)}
  where \eqn{K(r)} is the Ripley K-function of a spatially homogeneous
  point pattern, estimated by \code{\link{Kest}}.

  The inhomogeneous L-function is the corresponding transformation
  of the inhomogeneous K-function, estimated by \code{\link{Kinhom}}.
  It is appropriate when the point pattern clearly does not have a
  homogeneous intensity of points.

  The command \code{Linhom} first calls
  \code{\link{Kinhom}} to compute the estimate of the inhomogeneous K-function,
  and then applies the square root transformation.

  For a Poisson point pattern (homogeneous or inhomogeneous),
  the theoretical value of the inhomogeneous L-function is \eqn{L(r) = r}.
  The square root also has the effect of stabilising
  the variance of the estimator, so that L is more appropriate
  for use in simulation envelopes and hypothesis tests.
}

\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{L} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{L(r) = r}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{L(r)} obtained by the edge corrections
  named.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{Lest}},
  \code{\link{Kinhom}},
  \code{\link{pcf}}
}
\examples{
 data(japanesepines)
 X <- japanesepines
 L <- Linhom(X)
 plot(L)
 plot(L, main="L function for Japanese Pines")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{nonparametric}
