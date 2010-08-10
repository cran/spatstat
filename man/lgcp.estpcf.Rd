\name{lgcp.estpcf}
\alias{lgcp.estpcf}
\title{Fit a Log-Gaussian Cox Point Process by Minimum Contrast}
\description{
  Fits a log-Gaussian Cox point process model
  (with exponential covariance function)
  to a point pattern dataset by the Method of Minimum Contrast
  using the pair correlation function.
}
\usage{
lgcp.estpcf(X, startpar=c(sigma2=1,alpha=1), lambda=NULL,
q = 1/4, p = 2, rmin = NULL, rmax = NULL, ..., pcfargs=list())
}
\arguments{
  \item{X}{
    Data to which the model will be fitted.
    Either a point pattern or a summary statistic.
    See Details.
  }
  \item{startpar}{
    Vector of starting values for the parameters of the
    log-Gaussian Cox process model.
  }
  \item{lambda}{
    Optional. An estimate of the intensity of the point process.
  }
  \item{q,p}{
    Optional. Exponents for the contrast criterion.
  }
  \item{rmin, rmax}{
    Optional. The interval of \eqn{r} values for the contrast criterion.
  }
  \item{\dots}{
    Optional arguments passed to \code{\link[stats]{optim}}
    to control the optimisation algorithm. See Details.
  }
  \item{pcfargs}{
    Optional list containing arguments passed to \code{\link{pcf.ppp}}
    to control the smoothing in the estimation of the
    pair correlation function.
  }
}
\details{
  This algorithm fits a log-Gaussian Cox point process model
  to a point pattern dataset
  by the Method of Minimum Contrast, using the pair correlation function.

  The argument \code{X} can be either
  \describe{
    \item{a point pattern:}{An object of class \code{"ppp"}
      representing a point pattern dataset. 
      The pair correlation function of the point pattern will be computed
      using \code{\link{pcf}}, and the method of minimum contrast
      will be applied to this.
    }
    \item{a summary statistic:}{An object of class \code{"fv"} containing
      the values of a summary statistic, computed for a point pattern
      dataset. The summary statistic should be the pair correlation function,
      and this object should have been obtained by a call to
      \code{\link{pcf}} or one of its relatives.
    }
  }
  The algorithm fits a log-Gaussian Cox point process (LGCP)
  model to \code{X},  by finding the parameters of the LGCP model
  which give the closest match between the
  theoretical pair correlation function of the LGCP model
  and the observed pair correlation function.
  For a more detailed explanation of the Method of Minimum Contrast,
  see \code{\link{mincontrast}}.

  The model fitted is a stationary, isotropic log-Gaussian Cox process
  with exponential covariance
  (Moller and Waagepetersen, 2003, pp. 72-76).
  To define this process we start with
  a stationary Gaussian random field \eqn{Z} in the two-dimensional plane,
  with constant mean \eqn{\mu}{mu} and covariance function
  \deqn{
    c(r) = \sigma^2 e^{-r/\alpha}
  }{
    c(r) = sigma^2 * exp(-r/alpha)
  }
  where \eqn{\sigma^2}{sigma^2} and \eqn{\alpha}{alpha} are parameters.
  Given \eqn{Z}, we generate a Poisson point process \eqn{Y} with intensity
  function \eqn{\lambda(u) = \exp(Z(u))}{lambda(u) = exp(Z(u))} at
  location \eqn{u}. Then \eqn{Y} is a log-Gaussian Cox process.

  The theoretical pair correlation function of the LGCP is
  \deqn{
    g(r) = \exp(\sigma^2 \exp(-r/\alpha))
  }{
    g(r) = exp(sigma^2 * exp(-r/alpha))
  }
  The theoretical intensity of the LGCP is 
  \deqn{
    \lambda = \exp(\mu + \frac{\sigma^2}{2}).
  }{
    lambda= exp(mu + sigma^2/2).
  }
  In this algorithm, the Method of Minimum Contrast is first used to find
  optimal values of the parameters \eqn{\sigma^2}{sigma^2}
  and \eqn{\alpha}{alpha^2}. Then the remaining parameter
  \eqn{\mu}{mu} is inferred from the estimated intensity
  \eqn{\lambda}{lambda}.

  If the argument \code{lambda} is provided, then this is used
  as the value of \eqn{\lambda}{lambda}. Otherwise, if \code{X} is a
  point pattern, then  \eqn{\lambda}{lambda}
  will be estimated from \code{X}. 
  If \code{X} is a summary statistic and \code{lambda} is missing,
  then the intensity \eqn{\lambda}{lambda} cannot be estimated, and
  the parameter \eqn{\mu}{mu} will be returned as \code{NA}.

  The remaining arguments \code{rmin,rmax,q,p} control the
  method of minimum contrast; see \code{\link{mincontrast}}.

  The optimisation algorithm can be controlled through the
  additional arguments \code{"..."} which are passed to the
  optimisation function \code{\link[stats]{optim}}. For example,
  to constrain the parameter values to a certain range,
  use the argument \code{method="L-BFGS-B"} to select an optimisation
  algorithm that respects box constraints, and use the arguments
  \code{lower} and \code{upper} to specify (vectors of) minimum and
  maximum values for each parameter.
}
\value{
  An object of class \code{"minconfit"}. There are methods for printing
  and plotting this object. It contains the following main components:
  \item{par }{Vector of fitted parameter values.}
  \item{fit }{Function value table (object of class \code{"fv"})
    containing the observed values of the summary statistic
    (\code{observed}) and the theoretical values of the summary
    statistic computed from the fitted model parameters.
  }
}
\references{
  Moller, J. and Waagepetersen, R. (2003).
  Statistical Inference and Simulation for Spatial Point Processes.
  Chapman and Hall/CRC, Boca Raton.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\seealso{
  \code{\link{lgcp.estK}},
  \code{\link{thomas.estpcf}},
  \code{\link{matclust.estpcf}},
  \code{\link{mincontrast}},
  \code{\link{pcf}}
}
\examples{
    data(redwood)
    u <- lgcp.estpcf(redwood, c(sigma2=0.1, alpha=1))
    u
    plot(u)
}
\keyword{spatial}
\keyword{models}
