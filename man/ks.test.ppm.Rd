\name{ks.test.ppm}
\alias{ks.test.ppm}
\title{Kolmogorov-Smirnov Test for Point Process Model} 
\description{
  Performs a Kolmogorov-Smirnov test of goodness-of-fit
  of a Poisson point process model. The test compares the observed
  and predicted distributions of the values of a spatial covariate.
}
\usage{
ks.test.ppm(model, covariate, ...)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class \code{"ppm"}).
  }
  \item{covariate}{
    The spatial covariate on which the test will be based.
    An image (object of class \code{"im"})
    or a function.
  }
  \item{\dots}{
    Arguments passed to \code{\link{ks.test}}
    to control the test.
  }
}
\details{
  This function performs a goodness-of-fit test of
  a fitted point process model. The observed distribution
  of the values of a spatial covariate at the data points,
  and the predicted distribution of the same values under the model,
  are compared using the Kolmogorov-Smirnov test. 
  
  The argument \code{model} should be a fitted point process model
  (object of class \code{"ppm"}). It should be a Poisson point process.

  The argument \code{covariate} should be either a \code{function(x,y)}
  or a pixel image (object of class \code{"im"} containing the values
  of a spatial function.
  If \code{covariate} is an image, it should have numeric values,
  and its domain should cover the observation window of the
  \code{model}. If \code{covariate} is a function, it should expect
  two arguments \code{x} and \code{y} which are vectors of coordinates,
  and it should return a numeric vector of the same length
  as \code{x} and \code{y}.  

  First the original data point pattern is extracted from \code{model}.
  The values of the \code{covariate} at these data points are
  collected. 

  The predicted distribution of the values of the \code{covariate}
  under the fitted \code{model} is computed as follows.
  The values of the \code{covariate} at all locations in the
  observation window are evaluated,
  weighted according to the point process intensity of the fitted model,
  and compiled into a cumulative distribution function \eqn{F} using
  \code{\link{ewcdf}}.

  The probability integral transformation is then applied:
  the values of the \code{covariate} at the original data points
  are transformed by the predicted cumulative distribution function
  \eqn{F} into numbers between 0 and 1. If the model is correct,
  these numbers are i.i.d. uniform random numbers. The
  Kolmogorov-Smirnov test of uniformity is applied using
  \code{\link{ks.test}}.

  This test was apparently first described (in the context of
  spatial data) by Berman (1986). See also Baddeley et al (2005).
}
\value{
  An object of class \code{"htest"} containing the results of the
  test. See \code{\link{ks.test}} for details.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
  \code{\link{ks.test}},
  \code{\link{quadrat.test}},
  \code{\link{ppm}}
}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.

  Berman, M. (1986)
  Testing for spatial association between a point process
  and another stochastic process.
  \emph{Applied Statistics} \bold{35}, 54--62.
}
\examples{
   # nonuniform Poisson process
   X <- rpoispp(function(x,y) { 100 * exp(x) }, win=square(1))
   # fit uniform Poisson process
   fit0 <- ppm(X, ~1)
   # fit correct nonuniform Poisson process
   fit1 <- ppm(X, ~x)

   # test covariate = x coordinate
   xcoord <- function(x,y) { x }

   # test wrong model
   ks.test.ppm(fit0, xcoord)
   # test right model
   ks.test.ppm(fit1, xcoord)
}
\keyword{htest}
\keyword{spatial}

