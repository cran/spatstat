\name{suffstat}
\alias{suffstat}
\title{Sufficient Statistic of Point Process Model}
\description{
  The canonical sufficient statistic of a 
  point process model is evaluated for a given point pattern.
}
\usage{
  suffstat(model, X=data.ppm(model))
}
\arguments{
  \item{model}{A fitted point process model (object of class
    \code{"ppm"}).
  }
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
}
\value{
  A numeric vector of sufficient statistics. The entries
  correspond to the model coefficients \code{coef(model)}.
}
\details{
  The canonical sufficient statistic
  of \code{model} is evaluated for the point pattern \code{X}.
  This computation is useful for various Monte Carlo methods.
  
  Here \code{model} should be a point process model (object of class
  \code{"ppm"}, see \code{\link{ppm.object}}), typically obtained
  from the model-fitting function \code{\link{ppm}}. The argument
  \code{X} should be a point pattern (object of class \code{"ppp"}).

  Every point process model fitted by \code{\link{ppm}} has
  a probability density of the form
  \deqn{f(x) = Z(\theta) \exp(\theta^T S(x))}{f(x) = Z(theta) exp(theta * S(x))}
  where \eqn{x} denotes a typical realisation (i.e. a point pattern),
  \eqn{\theta}{theta} is the vector of model coefficients,
  \eqn{Z(\theta)}{Z(theta)} is a normalising constant,
  and \eqn{S(x)} is a function of the realisation \eqn{x}, called the
  ``canonical sufficient statistic'' of the model.

  For example, the stationary Poisson process has canonical sufficient
  statistic \eqn{S(x)=n(x)}, the number of points in \eqn{x}.
  The stationary Strauss process with interaction range \eqn{r}
  (and no edge correction) has canonical sufficient statistic
  \eqn{S(x)=(n(x),d(x))} where \eqn{d(x)} is the number of pairs
  of points in \eqn{x} which are closer than a distance \eqn{r}
  to each other. 

  \code{suffstat(model, X)} returns the value of \eqn{S(x)}, where \eqn{S} is
  the canonical sufficient statistic associated with \code{model},
  evaluated when \eqn{x} is the given point pattern \code{X}.
  The result is a numeric vector, with entries which correspond to the
  entries of the coefficient vector \code{coef(model)}.

  The sufficient statistic \eqn{S}
  does not depend on the fitted coefficients
  of the model. However it does depend on the irregular parameters
  which are fixed in the original call to \code{\link{ppm}}, for
  example, the interaction range \code{r} of the Strauss process.
  The sufficient statistic also depends on the edge correction that
  was used to fit the model.

  Non-finite values of the sufficient statistic (\code{NA} or
  \code{-Inf}) may be returned if the point pattern \code{X} is
  not a possible realisation of the model (i.e. if \code{X} has zero
  probability of occurring under \code{model} for all values of
  the canonical coefficients \eqn{\theta}{theta}).
}
\seealso{
  \code{\link{ppm}}
}
\examples{
    data(swedishpines)
    fitS <- ppm(swedishpines, ~1, Strauss(7))
    X <- rpoispp(summary(swedishpines)$intensity, win=swedishpines$window)
    suffstat(fitS, X)
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
