\name{MultiHard}
\alias{MultiHard}
\title{The Multitype Hard Core Point Process Model}
\description{
Creates an instance of the multitype hard core point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiHard(types, hradii)
}
\arguments{
  \item{types}{Vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
  \item{hradii}{Matrix of hard core radii}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype hard core process with
  hard core radii \eqn{hradii[i,j]}.
}
\details{
  This is a multitype version of the hard core process.
  A pair of points
  of types \eqn{i} and \eqn{j}
  must not lie closer than \eqn{h_{ij}}{h[i,j]} units apart.

  The matrix \code{hradii}
  must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no distance constraint should be applied
  for this combination of types.
  
  Note that only the hardcore radii
  are specified in \code{MultiHard}.
  The canonical parameters \eqn{\log(\beta_j)}{log(beta[j])}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{MultiHard()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{Strauss}}
}
\examples{
   h <- matrix(c(1,2,2,1), nrow=2,ncol=2)
   MultiHard(1:2, h)
   # prints a sensible description of itself
   data(betacells)
   h <- 15.0 * matrix(c(NA,1,1,NA), nrow=2,ncol=2)
   ppm(betacells, ~1, MultiHard(c("off","on"), h))
   # fit the stationary multitype hardcore process to `betacells'
}
\section{Warnings}{
  The argument \code{types} is interpreted as a
  set of factor levels. That is,
  in order that \code{\link{ppm}} can fit the multitype hard core model
  correctly to a point pattern \code{X},
  this must be a marked point pattern; the mark vector
  \code{marks(X)} must be a factor; and 
  the argument \code{types} must
  equal \code{levels(marks(X))}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
