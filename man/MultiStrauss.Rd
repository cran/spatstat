\name{MultiStrauss}
\alias{MultiStrauss}
\title{The Multitype Strauss Point Process Model}
\description{
Creates an instance of the multitype Strauss point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiStrauss(types, radii)
}
\arguments{
  \item{types}{Vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
  \item{radii}{Matrix of interaction radii}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype Strauss process with
  interaction radii \eqn{radii[i,j]}.
}
\details{
  The (stationary) multitype
  Strauss process with \eqn{m} types, with interaction radii
  \eqn{r_{ij}}{r[i,j]} and 
  parameters \eqn{\beta_j}{beta[j]} and \eqn{\gamma_{ij}}{gamma[i,j]}
  is the pairwise interaction point process
  in which each point of type \eqn{j}
  contributes a factor \eqn{\beta_j}{beta[j]} to the 
  probability density of the point pattern, and a pair of points
  of types \eqn{i} and \eqn{j} closer than \eqn{r_{ij}}{r[i,j]}
  units apart contributes a factor
  \eqn{\gamma_{ij}}{gamma[i,j]} to the density.

  The nonstationary multitype Strauss process is similar except that 
  the contribution of each individual point \eqn{x_i}{x[i]}
  is a function \eqn{\beta(x_i)}{beta(x[i])}
  of location and type, rather than a constant beta. 
 
  The function \code{\link{ppm}()}, which fits point process models to 
  point pattern data, requires an argument 
  of class \code{"interact"} describing the interpoint interaction
  structure of the model to be fitted. 
  The appropriate description of the multitype
  Strauss process pairwise interaction is
  yielded by the function \code{MultiStrauss()}. See the examples below.

  The matrix \code{radii} must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no interaction term should be included
  for this combination of types.
  
  Note that only the interaction radii are specified in \code{MultiStrauss}.
  The canonical parameters \eqn{\log(\beta_j)}{log(beta[j])}
  and \eqn{\log(\gamma_{ij})}{log(gamma[i,j])}
  are estimated by \code{\link{ppm}()}, not fixed in
  \code{Strauss()}. 
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{Strauss}}
}
\examples{
   r <- matrix(c(1,2,2,1), nrow=2,ncol=2)
   MultiStrauss(1:2, r)
   # prints a sensible description of itself
   data(ganglia)
   r <- 0.03 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   ppm(ganglia, ~1, MultiStrauss(c("off","on"), r), rbord=0.06)
   # fit the stationary multitype Strauss process to `ganglia'
   ppm(ganglia, ~polynom(x,y,3), MultiStrauss(c("off","on"), r), rbord=0.06)
   # fit a nonstationary Strauss process with log-cubic polynomial trend
}
\section{Warnings}{
  The argument \code{types} is interpreted as a
  set of factor levels. That is,
  in order that \code{\link{ppm}} can fit the multitype Strauss model
  correctly to a point pattern \code{X},
  this must be a marked point pattern; the mark vector
  \code{X$marks} must be a factor; and 
  the argument \code{types} must
  equal \code{levels(X$marks)}.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
