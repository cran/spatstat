\name{residuals.ppm}
\alias{residuals.ppm}
\title{
  Residuals for Fitted Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  compute residuals.
}
\synopsis{
  \method{residuals}{ppm}(object, type="raw", \dots, fittedvalues=fitted.ppm(object))
}
\usage{
  \method{residuals}{ppm}(object, type="raw", \dots)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which residuals should be calculated.
  }
  \item{type}{
    String indicating the type of residuals to be calculated.
    Current options are
    \code{"raw"}, \code{"inverse"} and \code{"pearson"}.
    A partial match is adequate.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A vector containing the (discretised) residuals.
  
  Entries in this vector correspond to the quadrature points (data or
  dummy points) used to fit the model. The quadrature points can be
  extracted from \code{object} by \code{union.quad(quad.ppm(object))}.
}
\details{
  This function computes several kinds of residuals for the fit of
  a point process model to a spatial point pattern dataset
  (Baddeley et al, 2005).
  Use \code{\link{diagnose.ppm}} to produce diagnostic plots
  based on these residuals.

  The argument \code{object} must be a fitted point process model
  (object of class \code{"ppm"}). Such objects are produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{ppm}}).
  This fitted model object contains complete
  information about the original data pattern. 

  Residuals are attached both to the data points and to some
  other points in the window of observation (namely, to the dummy
  points of the quadrature scheme used to fit the model).
  If the fitted model is correct, then the sum of the 
  residuals over all (data and dummy) points in a spatial region \eqn{B}
  has mean zero. For further explanation, see Baddeley et al (2005).
  
  The type of residual
  is chosen by the argument \code{type}. Current options are

  \describe{
    \item{\code{"raw"}:}{
    the raw residuals
    \deqn{
      r_j = z_j - w_j \lambda_j
    }{
      r[j] = z[j] - w[j] lambda[j]
    }
    at the quadrature points \eqn{u_j}{u[j]},
    where \eqn{z_j}{z[j]} is the indicator equal to 1 if \eqn{u_j}{u[j]}
    is a data point and 0 if \eqn{u_j}{u[j]} is a dummy point;
    \eqn{w_j}{w[j]} is the quadrature weight attached to
    \eqn{u_j}{u[j]}; and
    \deqn{\lambda_j = \hat\lambda(u_j,x)}{lambda[j] = lambda(u[j],x)}
    is the conditional intensity of the fitted model at \eqn{u_j}{u[j]}.
    These are the spatial analogue of the martingale residuals
    of a one-dimensional counting process.
  }
  \item{\code{"inverse"}:}{
    the `inverse-lambda' residuals (Baddeley et al, 2005)
    \deqn{
      r^{(I)}_j = \frac{r_j}{\lambda_j}
                   = \frac{z_j}{\lambda_j} - w_j
    }{
      rI[j] = r[j]/lambda[j] = z[j]/lambda[j] - w[j]
    }
    obtained by dividing the raw residuals by 
    the fitted conditional intensity. These are 
    a counterpart of the exponential energy marks (see \code{\link{eem}}).
  }
  \item{\code{"pearson"}:}{
    the Pearson residuals (Baddeley et al, 2005)
    \deqn{
      r^{(P)}_j = \frac{r_j}{\sqrt{\lambda_j}}
                   = \frac{z_j}{\sqrt{\lambda_j}}
                   - w_j \sqrt{\lambda_j}
    }{
      rP[j] = r[j]/sqrt(lambda[j])
      = z[j]/sqrt(lambda[j]) - w[j] sqrt(lambda[j])
    }
    obtained by dividing the raw residuals by the
    square root of the fitted conditional intensity.
    The Pearson residuals are standardised, in the sense
    that if the model (true and fitted) is Poisson,
    then the sum of the Pearson residuals in a spatial region \eqn{B}
    has variance equal to the area of \eqn{B}.
  }
  }

  Use \code{\link{diagnose.ppm}} to produce diagnostic plots
  based on these residuals.
}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\seealso{
 \code{\link{diagnose.ppm}},
 \code{\link{ppm.object}},
 \code{\link{ppm}}
}
\examples{
    data(cells)
    fit <- ppm(cells, ~x, Strauss(r=0.15), rbord=0.15)

    rp <- residuals.ppm(fit, type="pe")
    sum(rp) # should be close to 0

    # extract quadrature points to which the residuals correspond
    quadpoints <- union.quad(quad.ppm(fit))

    # plot residuals as marks attached to the quadrature points 
    quadmarked <- setmarks(quadpoints, rp)
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
\keyword{models}
\keyword{methods}
