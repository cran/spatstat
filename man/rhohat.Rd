\name{rhohat}
\alias{rhohat}
\title{
  Smoothing Estimate of Covariate Transformation
}
\description{
  Computes a smoothing estimate of the intensity of a point process,
  as a function of a spatial covariate.
}
\usage{
rhohat(object, covariate, ...,
       transform=FALSE,
       dimyx=NULL, eps=NULL,
       n = 512, bw = "nrd0", adjust=1, from = NULL, to = NULL,
       bwref=bw,
       covname)
}
\arguments{
  \item{object}{
    A point pattern (object of class \code{"ppp"}),
    a quadrature scheme (object of class \code{"quad"})
    or a fitted point process model (object of class \code{"ppm"}).
  }
  \item{covariate}{
    Either a \code{function(x,y)} or a pixel image (object of
    class \code{"im"}) providing the values of the covariate at any
    location.
    Alternatively one of the strings \code{"x"} or \code{"y"}
    signifying the Cartesian coordinates.
  }
  \item{transform}{
    Logical value determining the smoothing method. See Details.
  }
  \item{dimyx,eps}{
    Arguments passed to \code{\link{as.mask}} to control the pixel
    resolution at which the covariate will be evaluated.
  }
  \item{bw}{
    Smoothing bandwidth or bandwidth rule
    (passed to \code{\link{density.default}}).
  }
  \item{adjust}{
    Smoothing bandwidth adjustment factor
    (passed to \code{\link{density.default}}).
  }
  \item{n, from, to}{
    Arguments passed to \code{\link{density.default}} to
    control the number and range of values at which the function
    will be estimated.
  }
  \item{bwref}{
    Optional. An alternative value of \code{bw} to use when smoothing
    the reference density (the density of the covariate values
    observed at all locations in the window).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{density.default}}.
  }
  \item{covname}{
    Optional. Character string to use as the name of the covariate.
  }
}
\details{
  If \code{object} is a point pattern, this command assumes that
  \code{object} is a realisation of a Poisson point process with
  intensity function \eqn{\lambda(u)}{lambda(u)} of the form
  \deqn{\lambda(u) = \rho(Z(u))}{lambda(u) = rho(Z(u))} where
  \eqn{Z} is the spatial covariate function given by \code{covariate},
  and \eqn{\rho(z)}{rho(z)} is a function to be estimated. 
  This command computes an estimator of \eqn{\rho(z)}{rho(z)}
  proposed by Baddeley and Turner (2005).

  If \code{object} is a fitted point process model, suppose \code{X} is
  the original data point pattern to which the model was fitted. Then
  this command assumes \code{X} is a realisation of a Poisson point
  process with intensity function of the form
  \deqn{
    \lambda(u) = \rho(Z(u)) \kappa(u)
  }{
    lambda(u) = rho(Z(u)) * kappa(u)
  }
  where \eqn{\kappa(u)}{kappa(u)} is the intensity of the fitted model
  \code{object}. A modified version of the Baddeley-Turner (2005)
  smoothing estimator is computed.

  If \code{transform=FALSE}, the smoothing method is
  fixed bandwidth kernel smoothing, using \code{\link{density.default}}.
  If \code{transform=TRUE}, the smoothing method is variable-bandwidth kernel
  smoothing, implemented by applying the Probability Integral Transform
  to the covariate values, yielding values in the range 0 to 1,
  then applying edge-corrected fixed-bandwidth smoothing on the interval
  \eqn{[0,1]}, and back-transforming.
}
\value{
  A function value table (object of class \code{"fv"})
  containing the estimated values of \eqn{\rho}{rho} for a sequence
  of values of \eqn{Z}.
  Also belongs to the class \code{"rhohat"}
  which has special methods for \code{print}, \code{plot}
  and \code{predict}.
}
\references{
  Baddeley, A. and Chang, Y.-M. and Song, Y. and Turner, R. 
  \emph{Diagnostics for transformation of covariates
    in spatial Poisson point process models}.
  Submitted for publication.
  
  Baddeley, A. and Turner, R. (2005)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.  
}
\author{
  Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link[spatstat:methods.rhohat]{methods.rhohat}}
}
\examples{
  X <-  rpoispp(function(x,y){exp(3+3*x)})
  rho <- rhohat(X, "x")
  rho <- rhohat(X, function(x,y){x})
  plot(rho)
  curve(exp(3+3*x), lty=3, col=2, add=TRUE)
  fit <- ppm(X, ~x)
  rr <- rhohat(fit, "y")
}
\keyword{spatial}
\keyword{models}

