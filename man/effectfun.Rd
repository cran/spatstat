\name{effectfun}
\alias{effectfun}
\title{Compute Fitted Effect of a Spatial Covariate in a Point Process Model}
\description{
  Compute the intensity of a fitted point process model
  as a function of one of its covariates.
}
\usage{
  effectfun(model, covname, ...)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class
    \code{"ppm"}).
  }
  \item{covname}{
    The name of the covariate. A character string.
  }
  \item{\dots}{
    The fixed values of other covariates (in the form
    \code{name=value}) if required.
  }
}
\details{
  The object \code{model} should be an object of class
  \code{"ppm"} representing a point process model fitted to
  point pattern data.

  The model's trend formula should involve a spatial covariate
  named \code{covname}. This could be \code{"x"} or \code{"y"}
  representing one of the Cartesian coordinates.
  More commonly the covariate
  is another, external variable that was supplied when fitting the model.
    
  The command \code{effectfun} computes the fitted intensity
  of the point process \code{model} as a function of the covariate
  named \code{covname}.
  The return value can be plotted immediately, giving a
  plot of the fitted intensity against the value of the covariate.

  If the model also involves covariates other than \code{covname},
  then these covariates will be held fixed. Values for
  these other covariates must be provided as arguments
  to \code{effectfun} in the form \code{name=value}. 

  This command is just a wrapper for the prediction method
  \code{\link{predict.ppm}}. For more complicated computations
  about the fitted intensity, use \code{\link{predict.ppm}}. 
}
\value{
  A data frame containing a column of values of the covariate and a column
  of values of the fitted intensity.

  If the covariate named \code{covname} is numeric (rather than a factor
  or logical variable), the return value is
  also of class \code{"fv"} so that it can be plotted immediately.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{predict.ppm}},
  \code{\link{fv.object}}
}
\examples{
  data(copper)
  X <- copper$SouthPoints
  D <- distmap(copper$SouthLines)
  fit <- ppm(X, ~polynom(Z, 7), covariates=list(Z=D))
  plot(effectfun(fit, "Z"))
  fit <- ppm(X, ~x + polynom(Z, 7), covariates=list(Z=D))
  plot(effectfun(fit, "Z", x=20))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
