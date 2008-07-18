\name{kppm}
\alias{kppm}
\title{Fit cluster point process model}
\description{
  Fit a homogeneous or inhomogeneous cluster point process model
  to a point pattern.
}
\usage{
  kppm(X, trend = ~1, clusters = "Thomas", covariates = NULL, ...)
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}) to which the model
    should be fitted.
  }
  \item{trend}{
    An \R formula, with no left hand side,
    specifying the form of the log intensity.
  }
  \item{clusters}{
    Character string determining the cluster model.
    Partially matched.
    Options are \code{"Thomas"} and \code{"MatClust"}.
  }
  \item{covariates}{
    The values of any spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images.
  }
  \item{\dots}{
    Arguments passed to \code{\link{thomas.estK}} or
    \code{\link{matclust.estK}} controlling the minimum contrast
    fitting algorithm.
  }
}
\details{
  This function fits a cluster point process model to the
  point pattern dataset \code{X}.

  The algorithm first estimates the intensity function
  of the point process, by fitting a Poisson process with log intensity
  of the form specified by the formula \code{trend}.
  Then the inhomogeneous \eqn{K} function is estimated using the
  fitted intensity. Finally the parameters of the cluster model
  are estimated by the method of minimum contrast using the
  inhomogeneous \eqn{K} function.

  Currently the only options for the cluster mechanism
  are \code{clusters="Thomas"} for the Thomas process
  and \code{clusters="MatClust"} for the Matern cluster process.
}
\value{
  An object of class \code{"kppm"} representing the fitted model.
  There are methods for printing, plotting, predicting, simulating
  and updating objects of this class.
}
\seealso{
  \code{\link{plot.kppm}},
  \code{\link{predict.kppm}},
  \code{\link{simulate.kppm}},
  \code{\link{update.kppm}},
  \code{\link{thomas.estK}},
  \code{\link{matclust.estK}},
  \code{\link{mincontrast}},
  \code{\link{Kinhom}},
  \code{\link{ppm}}
}
\references{
  Waagepetersen, R. (2006). An estimation function approach to
  inference for inhomogeneous Neyman-Scott processes. Submitted.  
}
\examples{
  data(redwood)
  kppm(redwood, ~1, "Thomas")
  kppm(redwood, ~x, "MatClust") 
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

