\name{ppm.object}
\alias{ppm.object}
\title{Class of Fitted Point Process Models}
\description{
  A class \code{ppm} to represent a fitted stochastic model
  for a point process. The output of \code{\link{mpl}}.
}
\details{
  An object of class \code{ppm} represents a stochastic point process
  model that has been fitted to a point pattern dataset.
  Typically it is the output of the maximum pseudolikelihood model fitter,
  \code{\link{mpl}}. There are methods \code{\link{print.ppm}},
  \code{\link{predict.ppm}} and \code{\link{plot.ppm}}

  A \code{ppm} object contains at least the following entries:
  \tabular{ll}{
    \code{coef} \tab  the fitted regular parameters (as returned by
                     \code{glm}) \cr
    \code{trend} \tab  the trend formula or \code{NULL} \cr
    \code{interaction} \tab the point process interaction family 
                        (an object of class \code{"interact"})
			or \code{NULL} \cr
    \code{Q} \tab      the quadrature scheme used \cr
    \code{maxlogpl} \tab the maximised value of log pseudolikelihood \cr
    \code{correction} \tab  name of edge correction method used \cr
  }
  See \code{\link{mpl}} for explanation of these concepts.
  The irregular parameters (e.g. the interaction radius of the
  Strauss process) are encoded in the \code{interaction} entry.

  See also (for example) \code{\link{Strauss}}
  to understand how to specify 
  a point process model with unknown parameters.
}
\seealso{
  \code{\link{mpl}},
  \code{\link{print.ppm}},
  \code{\link{predict.ppm}},
  \code{\link{plot.ppm}}.
}
\section{Warnings}{
  The internal representation may change in the next few releases
  of the \code{spatstat} package.
}
\examples{
  \dontrun{library(spatstat)}
  data(cells)
  fit <- mpl(cells, ~ x, Strauss(0.1), correction="periodic")
  fit
  pred <- predict(fit)
  pred <- predict(fit, nx=50, ny=50, type="trend")
  \dontrun{
  plot(fit)
  }  
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
