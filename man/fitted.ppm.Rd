\name{fitted.ppm}
\alias{fitted.ppm}
\title{
  Fitted Conditional Intensity for Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  compute the fitted conditional intensity of the model
  at the points of the quadrature scheme used to fit the model.
}
\usage{
  fitted.ppm(object, \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A vector containing the values of the fitted conditional intensity.
  
  Entries in this vector correspond to the quadrature points (data or
  dummy points) used to fit the model. The quadrature points can be
  extracted from \code{object} by \code{union.quad(quad.ppm(object))}.
}
\details{
  The argument \code{object} must be a fitted point process model
  (object of class \code{"ppm"}). Such objects are produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{mpl}}).

  This function evaluates the conditional intensity
  \eqn{\hat\lambda(u, x)}{lambdahat(u,x)} of the fitted point process
  model for certain locations \eqn{u},
  where \code{x} is the original point pattern dataset to which
  the model was fitted.

  The locations \eqn{u} at which the fitted conditional intensity
  is evaluated, are the points of the
  quadrature scheme used to fit the model in \code{\link{mpl}}.
  They include the data points (the points of the original point pattern
  dataset \code{x}) and other ``dummy'' points 
  in the window of observation.

  Use \code{\link{predict.ppm}} to compute the fitted conditional
  intensity at other locations or with other values of the
  explanatory variables.
}
\references{
  Baddeley, A., Moller, J. and Turner, R. (2004)
  Residuals for spatial point processes.
  In preparation.
}
\seealso{
 \code{\link{ppm.object}},
 \code{\link{mpl}},
 \code{\link{predict.ppm}}
}
\examples{
    require(spatstat)
    data(cells)
    str <- mpl(cells, ~x, Strauss(r=0.15), rbord=0.15)
    lambda <- fitted(str)

    # extract quadrature points in corresponding order
    quadpoints <- union.quad(quad.ppm(str))

    # plot conditional intensity values
    # as circles centred on the quadrature points 
    quadmarked <- setmarks(quadpoints, lambda)
    plot(quadmarked)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
