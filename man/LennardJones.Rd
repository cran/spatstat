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
  Lennard-Jones pair potential (Lennard-Jones, 1924)
  each pair of points in the point pattern,
  a distance \eqn{d} apart,
  contributes a factor
  \deqn{
    \exp \left\{
    -
    4\epsilon
    \left[
       \left(
          \frac{\sigma}{d}
       \right)^{12}
       -
       \left(
          \frac{\sigma}{d}
       \right)^6
    \right]
    \right\}
  }{
    exp( - 4 * epsilon * ((sigma/d)^12 - (sigma/d)^6))
  }
  to the probability density,
  where \eqn{\sigma}{sigma} and \eqn{\epsilon}{epsilon} are
  positive parameters to be estimated.
  
  See \bold{Examples} for a plot of this expression.
  
  This potential causes very strong inhibition between points at short
  range, and attraction between points at medium range.
  The parameter  \eqn{\sigma}{sigma} is called the
  \emph{characteristic diameter} and controls the scale of interaction.
  The parameter \eqn{\epsilon}{epsilon} is called the \emph{well depth}
  and determines the strength of attraction.
  The potential switches from inhibition to attraction at
  \eqn{d=\sigma}{d=sigma}.
  The maximum value of the pair potential is
  \eqn{\exp(\epsilon)}{exp(epsilon)}
  occuring at distance
  \eqn{d = 2^{1/6} \sigma}{d = 2^(1/6) * sigma}.
  Interaction is usually considered to be negligible for distances
  \eqn{d > 2.5 \sigma \max\{1,\epsilon^{1/6}\}}{d > 2.5 * sigma * max(1, epsilon^(1/6))}.

  This potential is used 
  to model interactions between uncharged molecules in statistical physics.
  
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the Lennard-Jones pairwise interaction is
  yielded by the function \code{LennardJones()}.
  See the examples below.
 
  The ``canonical regular parameters'' estimated by \code{\link{ppm}} are
  \eqn{\theta_1 = 4 \epsilon \sigma^{12}}{theta1 = 4 * epsilon * sigma^12}
  and 
  \eqn{\theta_2 = 4 \epsilon \sigma^6}{theta2 = 4 * epsilon * sigma^6}.
}
\section{Warnings}{
  Fitting the Lennard-Jones model is extremely unstable, because
  of the strong dependence between the functions \eqn{d^{-12}}{d^(-12)}
  and \eqn{d^{-6}}{d^(-6)}. The fitting algorithm often fails to converge.
  
  The fitting algorithm may also fail to converge 
  when the spatial window of the point pattern dataset has dimensions
  that are much greater than 1 unit or much less than 1 unit,
  leading to overflow or underflow in the 12th power of the coordinates.
  To fix this problem, rescale the coordinates so that the window is
  close to the unit square (see Examples).
  
  Errors are likely to occur if this model is fitted to a point pattern dataset
  which does not exhibit both short-range inhibition and
  medium-range attraction between points.  The values of the parameters
  \eqn{\sigma}{sigma} and \eqn{\epsilon}{epsilon} may be \code{NA}
  (because the fitted canonical parameters have opposite sign, which
  usually occurs when the pattern is completely random).

  An absence of warnings does not mean that the fitted model is sensible.
  A negative value of \eqn{\epsilon}{epsilon} may be obtained (usually when
  the pattern is strongly clustered); this does not correspond
  to a valid point process model, but the software does not issue a warning.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\examples{
   data(demopat)
   demopat
   X <- rescale(unmark(demopat), 5000)
   X
   fit <- ppm(X, ~1, LennardJones(), rbord=0.1)
   fit
   plot(fitin(fit), xlim=c(0,0.01))
}
\references{
  Lennard-Jones, J.E. (1924) On the determination of molecular fields.
  \emph{Proc Royal Soc London A} \bold{106}, 463--477.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
