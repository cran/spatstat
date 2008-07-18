\name{plot.kppm}
\alias{plot.kppm}
\title{Plot a fitted cluster point process}
\description{
  Plots a fitted cluster point process model,
  displaying the fitted intensity and the fitted \eqn{K}-function.
}
\usage{
  \method{plot}{kppm}(x, ...)
}
\arguments{
  \item{x}{
    Fitted cluster point process model.
    An object of class \code{"kppm"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic function
  \code{\link{plot}} for the class \code{"kppm"} of fitted
  cluster point process models.
  
  The argument \code{x} should be a cluster point process model
  (object of class \code{"kppm"}) obtained using
  the function \code{\link{kppm}}.

  This command first plots the fitted intensity of the model,
  using \code{\link{plot.ppm}}, then displays the empirical and fitted
  inhomogeneous \eqn{K}-functions using \code{\link{plot.minconfit}}.
}
\value{
  Null.
}
\examples{
  data(redwood)
  fit <- kppm(redwood, ~1, "Thomas")
  plot(fit)
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
