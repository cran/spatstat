\name{MultiStraussHard}
\alias{MultiStraussHard}
\title{The Multitype/Hard Core Strauss Point Process Model}
\description{
Creates an instance of the multitype/hard core Strauss point process model
which can then be fitted to point pattern data.
}
\usage{
  MultiStraussHard(types, iradii, hradii)
}
\arguments{
  \item{types}{Vector of all possible types (i.e. the possible levels
    of the \code{marks} variable in the data)}
  \item{iradii}{Matrix of interaction radii}
  \item{hradii}{Matrix of hard core radii}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of the multitype/hard core Strauss process with
  interaction radii \eqn{iradii[i,j]}
  and hard core radii \eqn{hradii[i,j]}.
}
\details{
  This is a hybrid of the multitype Strauss process
  (see \code{\link{MultiStrauss}}) and the hard core process
  (case \eqn{\gamma=0}{gamma = 0} of the Strauss process).
  A pair of points
  of types \eqn{i} and \eqn{j}
  must not lie closer than \eqn{h_{ij}}{h[i,j]} units apart;
  if the pair lies more than \eqn{h_{ij}}{h[i,j]} and less than 
  \eqn{r_{ij}}{r[i,j]} units apart, it
  contributes a factor
  \eqn{\gamma_{ij}}{gamma[i,j]} to the probability density.

  The matrices \code{iradii} and \code{hradii}
  must be symmetric, with entries
  which are either positive numbers or \code{NA}. 
  A value of \code{NA} indicates that no interaction term should be included
  for this combination of types.
  
  Note that only the interaction radii and hardcore radii
  are specified in \code{MultiStraussHard}.
  The canonical parameters \eqn{\log(\beta_j)}{log(beta[j])}
  and \eqn{\log(\gamma_{ij})}{log(gamma[i,j])}
  are estimated by \code{\link{mpl}()}, not fixed in
  \code{MultiStraussHard()}. 
}
\seealso{
  \code{\link{mpl}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}},
  \code{\link{MultiStrauss}},
  \code{\link{Strauss}}
}
\examples{
   \dontrun{library(spatstat)}
   r <- matrix(3, nrow=2,ncol=2)
   h <- matrix(c(1,2,2,1), nrow=2,ncol=2)
   MultiStraussHard(1:2, r, h)
   # prints a sensible description of itself
   data(ganglia)
   r <- 0.03 * matrix(c(1,2,2,1), nrow=2,ncol=2)
   h <- 0.03 * matrix(c(NA,1,1,NA), nrow=2,ncol=2)
   mpl(ganglia, ~1, MultiStraussHard(c("off","on"), r, h), rbord=0.06)
   # fit the stationary multitype hardcore Strauss process to `ganglia'
}
\section{Warnings}{
  The argument \code{types} is interpreted as a
  set of factor levels. That is,
  in order that \code{\link{mpl}} can fit the multitype Strauss model
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
