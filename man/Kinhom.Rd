\name{Kinhom}
\alias{Kinhom}
\title{Inhomogeneous K-function}
\description{
  Estimates the inhomogeneous \eqn{K} function of
  a non-stationary point pattern.
}
\usage{
  Kinhom(X, lambda, r = NULL, breaks = NULL, slow = FALSE,
    correction=c("border", "isotropic", "Ripley", "translate"), \dots)
}
\arguments{
  \item{X}{
    The observed data point pattern,
    from which an estimate of the inhomogeneous \eqn{K} function
    will be computed.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
  }
  \item{lambda}{
    Vector of values of the estimated intensity function,
    evaluated at the points of the pattern \code{X}.
    Alternatively this may be a matrix: 
    see Details.
  }
  \item{r}{
    vector of values for the argument \eqn{r} at which
    the inhomogeneous \eqn{K} function
    should be evaluated. Not normally given by the user;
    there is a sensible default.
  }
  \item{breaks}{
    An alternative to the argument \code{r}.
    Not normally invoked by the user.
    See Details.
  }
  \item{slow}{
    Not normally given by the user. 
    Logical flag which selects the algorithm used to compute the
    inhomogeneous \eqn{K} function. The default (\code{slow=FALSE})
    is faster than the alternative (\code{slow=TRUE}). The slow
    algorithm is retained for internal purposes as a check on
    validity.
  }
  \item{correction}{
      A character vector containing any selection of the
      options \code{"border"}, \code{"bord.modif"},
      \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
      It specifies the edge correction(s) to be applied.
  }
  \item{\dots}{
    Currently ignored.
  }
} 
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing at least the following columns,
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{theo}{vector of values of \eqn{\pi r^2}{pi * r^2},
    the theoretical value of \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
    for an inhomogeneous Poisson process
  }
  and containing additional columns
  according to the choice specified in the \code{correction}
  argument. The additional columns are named
  \code{border}, \code{trans} and \code{iso}
  and give the estimated values of 
  \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
  using the border correction, translation correction,
  and Ripley isotropic correction, respectively.
}
\details{
  This computes a generalisation of the \eqn{K} function
  for inhomogeneous point patterns, proposed by
  Baddeley, Moller and Waagepetersen (2000).
  
  The ``ordinary'' \eqn{K} function
  (variously known as the reduced second order moment function
  and Ripley's \eqn{K} function), is
  described under \code{\link{Kest}}. It is defined only
  for stationary point processes.
  
  The inhomogeneous \eqn{K} function
  \eqn{K_{\rm inhom}(r)}{Kinhom(r)}
  is a direct generalisation to nonstationary point processes.
  Suppose \eqn{x} is a point process with non-constant intensity
  \eqn{\lambda(u)}{lambda(u)} at each location \eqn{u}.
  Define \eqn{K_{\rm inhom}(r)}{Kinhom(r)} to be the expected
  value, given that \eqn{u} is a point of \eqn{x},
  of the sum of all terms
  \eqn{1/\lambda(u)\lambda(x_j)}{1/lambda(u)lambda(x[j])}
  over all points \eqn{x_j}{x[j]}
  in the process separated from \eqn{u} by a distance less than \eqn{r}.
  This reduces to the ordinary \eqn{K} function if
  \eqn{\lambda()}{lambda()} is constant.
  If \eqn{x} is an inhomogeneous Poisson process with intensity
  function \eqn{\lambda(u)}{lambda(u)}, then
  \eqn{K_{\rm inhom}(r) = \pi r^2}{Kinhom(r) = pi * r^2}.

  This allows us to inspect a point pattern for evidence of 
  interpoint interactions after allowing for spatial inhomogeneity
  of the pattern. Values 
  \eqn{K_{\rm inhom}(r) > \pi r^2}{Kinhom(r) > pi * r^2}
  are suggestive of clustering.

  The argument \code{lambda} should be a vector of length equal to the
  number of points in the pattern \code{X}. It will be interpreted as
  giving the (estimated) values of \eqn{\lambda(x_i)}{lambda(x[i])} for
  each point \eqn{x_i}{x[i]} of the pattern \eqn{x}.

  Alternatively \code{lambda} may be a square matrix of dimensions
  \eqn{n \times n}{n * n} where \eqn{n} is the number of points in
  \code{X}. In this case \code{lambda[i,j]} will be interpreted
  as an estimate of
  \eqn{\lambda(x_i)\lambda(x_j)}{lambda(x[i]) * lambda(x[j])}.

  Edge corrections are used to correct bias in the estimation
  of \eqn{K_{\rm inhom}}{Kinhom}, in exactly the same way
  as for the classical \eqn{K} function.
  See the documentation for \code{\link{Kest}}.
  
  The pair correlation function can also be applied to the
  result of \code{Kinhom}; see \code{\link{pcf}}.
}
\references{
  Baddeley, A., Moller, J. and Waagepetersen, R. (2000)
  Non- and semiparametric estimation of interaction in
  inhomogeneous point patterns.
  \emph{Statistica Neerlandica} \bold{54}, 329--350.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{pcf}}
}
\examples{
  require(spatstat)

  data(lansing)
  # inhomogeneous pattern of maples
  X <- unmark(lansing[lansing$marks == "maple",])
  \testonly{
     sub <- sample(c(TRUE,FALSE), X$n, replace=TRUE, prob=c(0.1,0.9))
     X <- X[sub , ]
  }
  # fit spatial trend
  fit <- mpl(X, ~ polynom(x,y,2), Poisson())
  # predict intensity values at points themselves
  lambda <- predict(fit, newdata=data.frame(x=X$x, y=X$y), type="trend")
  # inhomogeneous K function
  Ki <- Kinhom(X, lambda)
  plot(Ki)

  # SIMULATED DATA
  # known intensity function
  lamfun <- function(x,y) { 100 * x }
  # inhomogeneous Poisson process
  Y <- rpoispp(lamfun, 100, owin())
  # evaluate intensity at points of pattern
  lambda <- lamfun(Y$x, Y$y)
  # inhomogeneous K function
  Ki <- Kinhom(Y, lambda)
  plot(Ki)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
