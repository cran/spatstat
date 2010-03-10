\name{logLik.ppm}
\alias{logLik.ppm}
\title{Log Likelihood for Poisson Point Process Model}
\description{
  Extracts the log likelihood of a fitted Poisson point process model.
}
\usage{
\method{logLik}{ppm}(object, ...)
}
\arguments{
  \item{object}{Fitted point process model.
    An object of class \code{"ppm"}.
  }
  \item{\dots}{Ignored.}
}
\details{
  The maximised value of the log likelihood for the fitted model
  (as approximated by quadrature using the Berman-Turner approximation)
  is extracted.

  If \code{object} is not a Poisson process, the maximised log
  \emph{pseudolikelihood} is returned, with a warning.
}
\value{
  A numerical value.
}
\seealso{
  \code{\link{ppm}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~x)
  logLik(fit)
  AIC(fit)
}
\keyword{spatial}
\keyword{models}

