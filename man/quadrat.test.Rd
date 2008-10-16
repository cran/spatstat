\name{quadrat.test}
\alias{quadrat.test}
\alias{quadrat.test.ppp}
\alias{quadrat.test.ppm}
\title{Chi-Squared Dispersion Test for Spatial Point Pattern Based on
  Quadrat Counts}
\description{
  Performs a chi-squared test of Complete Spatial Randomness
  for a given point pattern, based on quadrat counts.
  Alternatively performs a chi-squared goodness-of-fit test of a fitted
  inhomogeneous Poisson model.
}
\usage{
quadrat.test(X, ...)
\method{quadrat.test}{ppp}(X, nx=5, ny=nx, ..., xbreaks=NULL, ybreaks=NULL, tess=NULL)
\method{quadrat.test}{ppm}(X, nx=5, ny=nx, ..., xbreaks=NULL, ybreaks=NULL, tess=NULL)
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
  \item{\dots}{Ignored.}
  \item{xbreaks}{
    Optional. Numeric vector giving the \eqn{x} coordinates of the
    boundaries of the quadrats. Incompatible with \code{nx}.
  }
  \item{ybreaks}{
    Optional. Numeric vector giving the \eqn{y} coordinates of the
    boundaries of the quadrats. Incompatible with \code{ny}.
  }
  \item{tess}{
    Tessellation (object of class \code{"tess"}) determining the
    quadrats. Incompatible with \code{nx, ny, xbreaks, ybreaks}.
  }
}
\details{
  These functions perform \eqn{\chi^2}{chi^2} tests of goodness-of-fit
  for a point process model, based on quadrat counts.

  The function \code{quadrat.test} is generic, with methods for
  point patterns (class \code{"ppp"}) and point process models
  (class \code{"ppm"}).

  \itemize{
    \item
    if \code{X} is a point pattern, we test the null hypothesis
    that the data pattern is a realisation of Complete Spatial
    Randomness (the uniform Poisson point process).
    \item
    If \code{X} is a fitted point process model, then it should be
    a Poisson point process model. The 
    data to which this model was fitted are extracted from the model
    object, and are treated as the data point pattern for the test.
    We test the null hypothesis 
    that the data pattern is a realisation of the (inhomogeneous) Poisson point
    process specified by \code{X}.
  }

  In both cases, the window of observation is divided
  into tiles, and the number of data points in each tile is
  counted, as described in \code{\link{quadratcount}}.
  The quadrats are rectangular by default, or may be regions of arbitrary shape
  specified by the argument \code{tess}.
  The expected number of points in each quadrat is also calculated,
  as determined by CSR (in the first case) or by the fitted model
  (in the second case). Then we perform the
  \eqn{\chi^2}{chi^2} test of goodness-of-fit to the quadrat counts.

  The return value is an object of class \code{"htest"}.
  Printing the object gives comprehensible output
  about the outcome of the test.

  The return value also belongs to
  the special class \code{"quadrat.test"}. Plotting the object
  will display the quadrats, annotated by their observed and expected
  counts and the Pearson residuals. See the examples.
}
\seealso{
  \code{\link{quadratcount}},
  \code{\link{quadrats}},
  \code{\link{quadratresample}},
  \code{\link{chisq.test}},
  \code{\link{kstest}}.

  To test a Poisson point process model against a specific alternative,
  use \code{\link{anova.ppm}}.
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
  quadrat.test(simdat, 4, 3)

  # fitted model: inhomogeneous Poisson
  fitx <- ppm(simdat, ~x, Poisson())
  quadrat.test(fitx)

  te <- quadrat.test(simdat, 4)
  residuals(te)  # Pearson residuals

  plot(te)

  plot(simdat, pch="+", cols="green", cex=1.2, lwd=2)
  plot(te, add=TRUE, col="red", cex=1.5, lty=2, lwd=3)

  sublab <- eval(substitute(expression(p[chi^2]==z),
                       list(z=signif(te$p.value,3))))
  title(sub=sublab, cex.sub=3)

  # quadrats of irregular shape
  B <- dirichlet(runifpoint(6, simdat$window))
  qB <- quadrat.test(simdat, tess=B)
  plot(simdat, main="quadrat.test(simdat, tess=B)", pch="+")
  plot(qB, add=TRUE, col="red", lwd=2, cex=1.3)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
