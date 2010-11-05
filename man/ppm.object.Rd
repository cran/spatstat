\name{ppm.object}
\alias{ppm.object}
\title{Class of Fitted Point Process Models}
\description{
  A class \code{ppm} to represent a fitted stochastic model
  for a point process. The output of \code{\link{ppm}}.
}
\details{
  An object of class \code{ppm} represents a stochastic point process
  model that has been fitted to a point pattern dataset.
  Typically it is the output of the model fitter,
  \code{\link{ppm}}.

  There are methods \code{\link{print.ppm}},
  \code{\link{plot.ppm}}, \code{\link{predict.ppm}}, \code{\link{fitted.ppm}}
  and \code{\link{coef.ppm}} for the generic functions
  \code{\link{print}}, \code{\link{plot}}, \code{\link{predict}},
  \code{\link{fitted}} and \code{\link{coef}} respectively.

  See also (for example) \code{\link{Strauss}}
  to understand how to specify 
  a point process model with unknown parameters.

  Information about the data (to which the model was fitted)
  can be extracted using \code{\link{data.ppm}}, \code{\link{dummy.ppm}}
  and \code{\link{quad.ppm}}.

  If you really need to get at the internals,
  a \code{ppm} object contains at least the following entries:
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
  See \code{\link{ppm}} for explanation of these concepts.
  The irregular parameters (e.g. the interaction radius of the
  Strauss process) are encoded in the \code{interaction} entry.
  However see the Warnings.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{coef.ppm}},
  \code{\link{fitted.ppm}},
  \code{\link{print.ppm}},
  \code{\link{predict.ppm}},
  \code{\link{plot.ppm}}.
}
\section{Warnings}{
  The internal representation may change in the next few releases
  of the \pkg{spatstat} package.
}
\examples{
  data(cells)
  fit <- ppm(cells, ~ x, Strauss(0.1), correction="periodic")
  fit
  coef(fit)
  \dontrun{
  pred <- predict(fit)
  }
  pred <- predict(fit, ngrid=20, type="trend")
  \dontrun{
  plot(fit)
  }  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{attribute}
 
 
