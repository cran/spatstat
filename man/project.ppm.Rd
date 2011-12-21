\name{project.ppm}
\alias{project.ppm}
\title{
  Force Point Process Model to be Valid
}
\description{
  Ensures that a fitted point process model 
  satisfies the integrability conditions for existence of the point process.
}
\usage{
project.ppm(object, fatal=FALSE)
}
\arguments{
  \item{object}{
    Fitted point process model (object of class \code{"ppm"}).
  }
  \item{fatal}{
    Logical. Whether to generate an error if the model cannot be
    projected to a valid model.
  }
}
\details{
  The model-fitting function \code{\link{ppm}}
  fits Gibbs point process models to point pattern data.
  By default, the fitted model returned by \code{\link{ppm}} may not
  actually exist as a point process. 

  Unlike a regression model, which is well-defined for any values
  of the fitted regression coefficients, a Gibbs point process model
  is only well-defined if the fitted interaction parameters 
  satisfy some constraints. 
  A famous example is the Strauss process (see \code{\link{Strauss}})
  which exists only when the interaction parameter \eqn{\gamma}{gamma}
  is less than or equal to 1. For values \eqn{\gamma > 1}{gamma > 1},
  the probability density is not integrable and the process does not
  exist (and cannot be simulated).

  By default, \code{\link{ppm}} does not enforce the constraint that
  a fitted Strauss process (for example) must satisfy
  \eqn{\gamma \le 1}{gamma <= 1}.
  This is because a fitted parameter value of \eqn{\gamma > 1}{gamma > 1} 
  could be useful information for data analysis, as it indicates that
  the Strauss model is not appropriate, and suggests a clustered model should be
  fitted.

  The function \code{project.ppm} changes the coefficients
  in the fitted model \code{object} so that the model is valid.
  Because the loglikelihood and log pseudolikelihood are concave,
  this simply involves projecting the coefficients onto the boundary
  of the parameter space.

  Use the function \code{\link{valid.ppm}} to check whether a fitted model
  object specifies a well-defined point process.
}
\value{
  Another point process model (object of class \code{"ppm"}).
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppm}},
  \code{\link{project.ppm}}
}
\examples{
   fit <- ppm(redwood, ~1, Strauss(0.1))
   coef(fit)
   fit2 <- project.ppm(fit)
   coef(fit2)
}
\keyword{spatial}
\keyword{models}
