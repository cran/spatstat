\name{vcov.ppm}
\alias{vcov.ppm}
\title{Variance-Covariance Matrix for a Fitted Point Process Model}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a fitted point process model.
}
\usage{
  \method{vcov}{ppm}(object, \dots, what = "vcov", verbose = TRUE,
gam.action="warn", matrix.action="warn", hessian=FALSE)
}
\arguments{
  \item{object}{A fitted point process model (an object of class \code{"ppm"}.)}
  \item{\dots}{Ignored.}
  \item{what}{Character string (partially-matched)
    that specifies what matrix is returned.
    Options are \code{"vcov"} for the variance-covariance matrix,
    \code{"corr"} for the correlation matrix, and
    \code{"fisher"} or \code{"Fisher"}
    for the Fisher information matrix.
  }
  \item{verbose}{Logical. If \code{TRUE}, a message will be printed
    if various minor problems are encountered.
  }
  \item{gam.action}{String indicating what to do if \code{object} was
    fitted by \code{gam}. Options are \code{"fatal"}, \code{"warn"} and
    \code{"silent"}.
  }
  \item{matrix.action}{String indicating what to do if the matrix
    is ill-conditioned (so that its inverse cannot be calculated).
    Options are \code{"fatal"}, \code{"warn"} and
    \code{"silent"}.
  }
  \item{hessian}{
    Logical. Use the negative Hessian matrix
    of the log pseudolikelihood instead of the Fisher information.
  }
}
\details{
  This function computes the asymptotic variance-covariance
  matrix of the estimates of the canonical parameters in the
  point process model \code{object}. It is a method for the 
  generic function \code{\link{vcov}}.

  \code{object} should be an object of class \code{"ppm"}, typically
  produced by \code{\link{ppm}}.

  The canonical parameters of the fitted model \code{object}
  are the quantities returned by \code{coef.ppm(object)}.
  The function \code{vcov} calculates the variance-covariance matrix
  for these parameters.
  
  The argument \code{what} provides three options:
  \describe{
    \item{\code{what="vcov"}}{
      return the variance-covariance matrix of the parameter estimates
    }
    \item{\code{what="corr"}}{
      return the correlation matrix of the parameter estimates
    }
    \item{\code{what="fisher"}}{
      return the observed Fisher information matrix.
    }
  }
  In all three cases, the result is a square matrix.
  The rows and columns of the matrix correspond to the canonical
  parameters given by \code{\link{coef.ppm}(object)}. The row and column
  names of the matrix are also identical to the names in
  \code{\link{coef.ppm}(object)}.

  For models fitted by maximum pseudolikelihood (which is the
  default in \code{\link{ppm}}), the current implementation only works
  for Poisson point processes. 
  The calculations are based on standard asymptotic theory for the maximum
  likelihood estimator.
  The observed Fisher information matrix of the fitted model
  \code{object} is first computed, by
  summing over the Berman-Turner quadrature points in the fitted model.
  The asymptotic variance-covariance matrix is calculated as the inverse of the
  observed Fisher information. The correlation matrix is then obtained
  by normalising.

  For models fitted by the Huang-Ogata method (\code{method="ho"} in
  the call to \code{\link{ppm}}), the implementation works for all
  models. A Monte Carlo estimate of the Fisher information matrix is
  calculated using the results of the original fit. 
  
  The argument \code{verbose} makes it possible to suppress some
  diagnostic messages.

  The asymptotic theory is not correct if the model was fitted using
  \code{gam} (by calling \code{\link{ppm}} with \code{use.gam=TRUE}).
  The argument \code{gamaction} determines what to do in this case.
  If \code{gamaction="fatal"}, an error is generated.
  If \code{gamaction="warn"}, a warning is issued and the calculation
  proceeds using the incorrect theory for the parametric case, which is
  probably a reasonable approximation in many applications.
  If \code{gamaction="silent"}, the calculation proceeds without a
  warning.
  
  If \code{hessian=TRUE} then the negative Hessian (second derivative)
  matrix of the log pseudolikelihood, and its inverse,
  will be computed.
  For non-Poisson models, this is not a valid estimate of variance,
  but is useful for other calculations.
}
\section{Error messages}{
  An error message that reports
  \emph{system is computationally singular} indicates that the
  determinant of the Fisher information matrix was either too large 
  or too small for reliable numerical calculation.
  This can occur either because of numerical overflow or
  because of collinearity in the covariates. Most commonly it occurs
  because of numerical overflow: to check this, rescale the 
  coordinates of the data points and refit the model. See the Examples.
}
\value{
  A square matrix.
}
\examples{
  X <- rpoispp(42)
  fit <- ppm(X, ~ x + y)
  vcov(fit)
  vcov(fit, what="Fish")

  # example of singular system
  data(demopat)
  m <- ppm(demopat, ~polynom(x,y,2))
  \dontrun{
    try(v <- vcov(m))
  }
  # rescale x, y coordinates to range [0,1] x [0,1] approximately
  demopat <- rescale(demopat, 10000)
  m <- ppm(demopat, ~polynom(x,y,2))
  v <- vcov(m)  
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{models}

