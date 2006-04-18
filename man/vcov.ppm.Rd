\name{vcov.ppm}
\alias{vcov.ppm}
\title{Variance-Covariance Matrix for a Fitted Point Process Model}
\description{
  Returns the variance-covariance matrix of the estimates of the
  parameters of a fitted point process model.
}
\usage{
  \method{vcov}{ppm}(object, \dots, what = "vcov", verbose = TRUE)
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
    if various minor problems are encountered.}
  
}
\details{
  This function computes the asymptotic variance-covariance
  matrix of the estimates of the canonical parameters in the
  point process model \code{object}. It is a method for the 
  generic function \code{\link{vcov}}.

  \code{object} should be an object of class \code{"ppm"}, typically
  produced by \code{\link{ppm}}. The current implementation only works
  for Poisson point processes.

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
  The calculations are based on standard asymptotic theory for the maximum
  likelihood estimator.
  In all cases, the observed Fisher information matrix of the fitted model
  \code{object} is first computed, by
  summing over the Berman-Turner quadrature points in the fitted model.
  The asymptotic variance-covariance matrix is calculated as the inverse of the
  observed Fisher information. The correlation matrix is then obtained
  by normalising.

  In all three cases, the result is a square matrix.
  The rows and columns of the matrix correspond to the canonical
  parameters given by \code{\link{coef.ppm}(object)}. The row and column
  names of the matrix are also identical to the names in
  \code{\link{coef.ppm}(object)}.

  The argument \code{verbose} makes it possible to suppress some
  diagnostic messages.
}
\value{
  A square matrix.
}
\examples{
  X <- rpoispp(42)
  fit <- ppm(X, ~ x + y)
  vcov(fit)
  vcov(fit, what="Fish")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
