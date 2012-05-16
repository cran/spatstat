\name{simulate.ppm}
\alias{simulate.ppm}
\title{Simulate a Fitted Gibbs Point Process Model}
\description{
  Generates simulated realisations from a fitted Gibbs or Poisson
  point process model.
}
\usage{
  \method{simulate}{ppm}(object, nsim=1, ...,
                         start = NULL,
                         control = default.rmhcontrol(object),
                         project=TRUE,
                         verbose=FALSE, progress=(nsim > 1))
}
\arguments{
  \item{object}{
    Fitted point process model.
    An object of class \code{"ppm"}.
  }
  \item{nsim}{
    Number of simulated realisations.
  }
  \item{start}{Data determining the initial state
    of the Metropolis-Hastings algorithm.  See
    \code{\link{rmhstart}} for description of these arguments.
    Defaults to \code{list(x.start=data.ppm(model))}
  }
  \item{control}{Data controlling the running of
    the Metropolis-Hastings algorithm.  See \code{\link{rmhcontrol}}
    for description of these arguments.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{rmhcontrol}},
    or to \code{\link{rmh.default}}, or to covariate functions in the model.
  }
  \item{project}{
    Logical flag indicating what to do if the fitted model is
    invalid (in the sense that the values of the fitted coefficients do not
    specify a valid point process).
    If \code{project=TRUE} the closest valid model will be simulated;
    if \code{project=FALSE} an error will occur.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports
    from \code{\link{rmh.ppm}}
    during the simulation of each point pattern.
  }
  \item{progress}{
    Logical flag indicating whether to print progress reports
    for the sequence of simulations.
  }
}
\details{
  This function is a method for the generic function
  \code{\link[stats]{simulate}} for the class \code{"ppm"} of fitted
  point process models.
  
  Simulations are performed by \code{\link{rmh.ppm}}.
}
\value{
  A list of length \code{nsim} containing simulated point patterns
  (objects of class \code{"ppp"}).
}
\examples{
  \testonly{op <- spatstat.options(rmh.nrep=10)}
  fit <- ppm(japanesepines, ~1, Strauss(0.1))
  simulate(fit, 2)
  \testonly{spatstat.options(op)}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{simulate.kppm}},
  \code{\link[stats]{simulate}}
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
