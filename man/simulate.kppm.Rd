\name{simulate.kppm}
\alias{simulate.kppm}
\title{Simulate a fitted cluster point process model.}
\description{
  Generates simulated realisations from a fitted cluster point process model.
}
\usage{
  \method{simulate}{kppm}(object, nsim = 1, ...)
}
\arguments{
  \item{object}{
    Fitted cluster point process model. An object of class \code{"kppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{\dots}{Ignored.}
}
\details{
  This function is a method for the generic function
  \code{\link{simulate}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  Simulations are performed by \code{\link{rThomas}}
  or \code{\link{rMatClust}} depending on the cluster mechanism.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}).
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~1, "Thomas")
  simulate(fit, 2)
}
\seealso{
  \code{\link{kppm}},
  \code{\link{plot.ppm}},
  \code{\link{plot.minconfit}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
