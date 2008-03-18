\name{Lest}
\alias{Lest}
\title{L-function}
\description{
  Calculates an estimate of Ripley's L-function
  for a spatial point pattern.
}
\usage{
  Lest(...)
}
\arguments{
  \item{\dots}{
    Arguments passed to \code{\link{Kest}}
    to estimate the K-function.
  }
}
\details{
  This command computes an estimate of the L-function
  for a spatial point pattern.
  The L-function is a transformation of Ripley's K-function,
  \deqn{L(r) = \sqrt{\frac{K(r)}{\pi}}}{L(r) = sqrt(K(r)/pi)}
  where \eqn{K(r)} is the K-function.

  See \code{\link{Kest}} for information
  about Ripley's K-function.

  The command \code{Lest} first calls
  \code{\link{Kest}} to compute the estimate of the K-function,
  and then applies the square root transformation.

  For a completely random (uniform Poisson) point pattern,
  the theoretical value of the L-function is \eqn{L(r) = r}.
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
  \code{\link{pcf}}
}
\examples{
 data(cells)
 L <- Lest(cells)
 plot(L)
 plot(L, main="L function for cells")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
