\name{LennardJones}
\alias{LennardJones}
\title{The Lennard-Jones Potential}
\description{
   Creates the Lennard-Jones pairwise interaction structure
   which can then be fitted to point pattern data.
}
\usage{
  LennardJones()
}
\value{
  An object of class \code{"interact"}
  describing the Lennard-Jones interpoint interaction
  structure.
}
\details{
  In a pairwise interaction point process with the
  Lennard-Jones pair potential, 
  each pair of points in the point pattern,
  a distance \eqn{d} apart,
  contributes a factor
  \deqn{
    \exp \left\{
	- 
          \left(
            \frac{\sigma}{d}
          \right)^{12}
    + \tau
          \left(
          \frac{\sigma}{d}
          \right)^6
    \right\}
  }{
    exp( - (sigma/d)^12 + tau (sigma/d)^6)
  }
  to the probability density,
  where \eqn{\sigma}{sigma} and \eqn{\tau}{tau} are
  positive parameters to be estimated. 

  See \bold{Examples} for a plot of this expression.
  
  This potential causes very strong inhibition between points at short range,
  and attraction between points at medium range.
  Roughly speaking, \eqn{\sigma}{sigma} controls the scale of both
  types of interaction, and \eqn{\tau}{tau} determines the strength of 
  attraction.
  The potential switches from inhibition to attraction at
  \eqn{d=\sigma/\tau^{1/6}}{d=sigma/tau^(1/6)}.
  Maximum attraction occurs at distance
  \eqn{d = (\frac 2 \tau)^{1/6} \sigma}{d = (2/tau)^(1/6) sigma}
  and the maximum achieved is \eqn{\exp(\tau^2/4)}{exp((tau^2)/4)}.
  Interaction is negligible for distances
  \eqn{d > 2 \sigma \max\{1,\tau^{1/6}\}}{d > 2 * sigma * max(1, tau^(1/6))}.

  This potential 
  is used (in a slightly different parameterisation)
  to model interactions between uncharged molecules in statistical physics.
  
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Lennard-Jones pairwise interaction is
  yielded by the function \code{LennardJones()}.
  See the examples below.
 
  The ``canonical regular parameters'' estimated by \code{\link{ppm}} are
  \eqn{\theta_1 = \sigma^{12}}{theta1 = sigma^12}
  and 
  \eqn{\theta_2 = \tau \sigma^6}{theta2 = tau * sigma^6}.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\examples{
   X <- rpoispp(100)
   ppm(X, ~1, LennardJones(), correction="translate")
   # Typically yields very small values for theta_1, theta_2
   # so the values of sigma, tau may not be sensible
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
