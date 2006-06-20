\name{quadrat.test}
\alias{quadrat.test}
\title{Chi-Squared Dispersion Test for Spatial Point Pattern Based on
  Quadrat Counts}
\description{
  Performs a chi-squared test of complete spatial randomness
  for a given point pattern, based on quadrat counts.
  Alternatively performs a chi-squared goodness-of-fit test of a fitted
  inhomogeneous Poisson model.
}
\usage{
quadrat.test(X, nx = 5, ny = nx, xbreaks = NULL, ybreaks = NULL, fit)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"})
    to be subjected to the goodness-of-fit test.
    Alternatively a fitted point process model (object of class
    \code{"ppm"}) to be tested.
    }
    \item{nx,ny}{
      Numbers of quadrats in the \eqn{x} and \eqn{y} directions.
      Incompatible with \code{xbreaks} and \code{ybreaks}.
    }
    \item{xbreaks}{
      Optional. Numeric vector giving the \eqn{x} coordinates of the
      boundaries of the quadrats. Incompatible with \code{nx}.
    }
    \item{ybreaks}{
      Optional. Numeric vector giving the \eqn{y} coordinates of the
      boundaries of the quadrats. Incompatible with \code{ny}.
    }
    \item{fit}{
      Optional. A fitted point process model (object of class
      \code{"ppm"}). The point pattern \code{X} will be subjected to
      a test of goodness-of-fit to the model \code{fit}.
    }
}
\details{
  This function performs a \eqn{\chi^2}{chi^2} test of goodness-of-fit
  to the Poisson point process (including \sQuote{Complete Spatial
    Randomness} but also inhomogeneous Poisson processes),
  based on quadrat counts.

  If \code{X} is a point pattern, it is taken as the data point pattern
  for the test. If \code{X} is a fitted point process model, then the
  data to which this model was fitted are extracted from the model
  object, and are treated as the data point pattern for the test.
  
  The window of observation is divided into rectangular tiles
  and the number of data points in each tile is counted,
  as described in \code{\link{quadratcount}}.

  If \code{fit} is absent, then we test the null hypothesis
  that the data pattern is a realisation of Complete Spatial
  Randomness (the uniform Poisson point process) by applying the
  \eqn{\chi^2}{chi^2} test of goodness-of-fit to the quadrat counts.

  If \code{fit} is present, then it should be a point process model
  (object of class \code{"ppm"}) and it should be a Poisson point
  process. Then we test the null hypothesis 
  that the data pattern is a realisation of the (inhomogeneous) Poisson point
  process specified by \code{fit}. Again this is a 
  \eqn{\chi^2}{chi^2} test of goodness-of-fit to the quadrat counts.

  The return value is an object of class \code{"htest"}.
  Printing the object gives comprehensible output
  about the outcome of the test. The return value also belongs to
  the special class \code{"quadrat.test"}. Plotting the object
  will display the quadrats, annotated by their observed and expected
  counts and the Pearson residuals. See the examples.

  To test the Poisson point process against a specific alternative
  point process model, use \code{\link{anova.ppm}}.
}
\seealso{
  \code{\link{quadratcount}},
  \code{\link{chisq.test}}
}
\value{
  An object of class \code{"htest"}. See \code{\link{chisq.test}}
  for explanation.

  The return value is also an object of the special class
  \code{"quadrat.test"}, and there is a plot method for this class.
  See the examples.
}
\examples{
  data(simdat)
  quadrat.test(simdat)
  quadrat.test(simdat, 4)

  # fitted model: inhomogeneous Poisson
  fitx <- ppm(simdat, ~x, Poisson())
  # equivalent:
  quadrat.test(simdat, fit=fitx)
  quadrat.test(fitx)

  te <- quadrat.test(simdat, 4)
  residuals(te)  # Pearson residuals

  plot(te)

  plot(simdat, pch="+", col="green", cex=1.2, lwd=2)
  plot(te, add=TRUE, col="red", cex=1.5, lty=2, lwd=3)

  sublab <- eval(substitute(expression(p[chi^2]==z),
                       list(z=signif(te$p.value,3))))
  title(sub=sublab, cex.sub=3)

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{htest}
