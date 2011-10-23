\name{default.expand}
\alias{default.expand}
\title{Compute Expansion Window for Simulation}
\description{
  Defines the default expansion window
  for simulation of a fitted point process model.
}
\usage{
  default.expand(object, m=2, epsilon=1e-6)
}
\arguments{
  \item{object}{
    A fitted point process model (object of class \code{"ppp"}).
  }
  \item{m}{
    A single numeric value.
    The window will be expanded by a distance 
    \code{m * reach(object)} along each side.
  }
  \item{epsilon}{
    Threshold argument passed to \code{\link{reach.ppm}} to determine
    \code{reach(object)}.
  }
} 
\value{
  A window (object of class \code{"owin"})
  or the numeric value 1, representing no expansion.
}
\details{
  This function computes a default value for the
  expansion window (the argument \code{expand} in \code{\link{rmhcontrol}})
  given a fitted point process model \code{object}.
  This default is used by \code{\link{envelope}}, 
  \code{\link{qqplot.ppm}} and other functions.

  Suppose we wish to generate simulated realisations
  of a fitted point process model inside a window \code{w}.
  It is advisable to first simulate
  the pattern on a larger window, and then clip it to the original
  window \code{w}. This avoids edge effects in the simulation.
  It is called \emph{expansion} of the simulation window.
  
  Accordingly, for the Metropolis-Hastings simulation algorithm
  \code{\link{rmh}}, the algorithm control parameters specified by
  \code{\link{rmhcontrol}} include an argument \code{expand} that
  determines the expansion of the simulation window.

  The function \code{default.expand} determines the default expansion
  window for a fitted point process model \code{object}.
  
  If the model is Poisson, then no expansion is necessary.
  No expansion is performed by default,
  and \code{default.expand} returns the value 1, representing no expansion.
  The simulation window is the original window \code{w = as.owin(object)}.
  
  If the model depends on external covariates (i.e.\ covariates other than
  the Cartesian covariates \code{x} and \code{y} and the \code{marks})
  then no expansion is feasible, in general, because the spatial domain
  of the covariates is not guaranteed to be large enough.
  \code{default.expand} returns the value 1, representing no expansion.
  The simulation window is the original window \code{w = as.owin(object)}.

  Otherwise, the original window \code{w = as.owin(object)} is expanded by
  a distance \code{m * rr}, where
  \code{rr} is the interaction range of the model, computed by
  \code{\link{reach}}. If \code{w} is a rectangle then
  each edge of \code{w} is displaced outward by distance \code{m * rr}.
  If \code{w} is not a rectangle then \code{w} is dilated by
  distance \code{m * rr} using \code{\link{dilation}}.
}
\seealso{
  \code{\link{rmhcontrol}},
  \code{\link{rmh}},
  \code{\link{envelope}},
  \code{\link{qqplot.ppm}}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~1, Strauss(0.07))
  default.expand(fit)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
