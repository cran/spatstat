\name{plot.kstest}
\alias{plot.kstest}
\title{Plot a Spatial Kolmogorov-Smirnov Test}
\description{
  Plot the result of a spatial Kolmogorov-Smirnov test
}
\usage{
\method{plot}{kstest}(x, ...,
                   lwd=par("lwd"), col=par("col"), lty=par("lty"),
                   lwd0=lwd, col0=col, lty0=lty)
}
\arguments{
  \item{x}{
    Object to be plotted. An object of class \code{"kstest"}
    produced by a method for \code{\link{kstest}}.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting function
    \code{\link{plot.default}}.
  }
  \item{col,lwd,lty}{
    The width, colour and type of lines used to plot the
    empirical distribution.
  }
  \item{col0,lwd0,lty0}{
    The width, colour and type of lines used to plot the
    predicted distribution.
  }
}
\value{
  \code{NULL}.
}
\details{
  This is the \code{plot} method for the class \code{"kstest"}.
  An object of this class represents the outcome of
  a spatial Kolmogorov-Smirnov test, computed by \code{\link{kstest}}.

  The plot displays the two cumulative distribution functions
  that are compared by the test: namely the empirical cumulative distribution
  function of the covariate at the data points, and the predicted
  cumulative distribution function of the covariate under the model,
  both plotted against the value of the covariate.
}
\seealso{
  \code{\link{kstest}}
}
\examples{
   # synthetic data: nonuniform Poisson process
   X <- rpoispp(function(x,y) { 100 * exp(x) }, win=square(1))

   # fit uniform Poisson process
   fit0 <- ppm(X, ~1)

   # test covariate = x coordinate
   xcoord <- function(x,y) { x }

   # test wrong model
   k <- kstest(fit0, xcoord)

   # plot result of test
   plot(k, lwd0=3)
}

\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
