\name{envelope.envelope}
\alias{envelope.envelope}
\title{
  Recompute Envelopes
}
\description{
  Given a simulation envelope (object of class \code{"envelope"}),
  compute another envelope from the same simulation data
  using different parameters.
}
\usage{
\method{envelope}{envelope}(Y, fun = NULL, ...)
}
\arguments{
  \item{Y}{
    A simulation envelope (object of class \code{"envelope"}).
}
  \item{fun}{
    Optional. Summary function to be applied to the simulated point patterns.
  }
  \item{\dots}{
    Parameters controlling the type of envelope that is re-computed. 
    See \code{\link{envelope}}.
  }
}
\details{
  This function can be used to re-compute a simulation envelope
  from previously simulated data, using different parameter settings
  for the envelope: for example, a different
  significance level, or a global envelope instead of a pointwise
  envelope.

  The function \code{\link{envelope}} is generic. This is the method for
  the class \code{"envelope"}. 
  
  The argument \code{Y} should be a simulation envelope (object of
  class \code{"envelope"}) produced by any of the methods for
  \code{\link{envelope}}.

  If the argument \code{fun} is given, it should be a summary function
  that can be applied to the simulated point patterns that were
  used to create \code{Y}. Also, \code{Y} should contain
  these point patterns, so \code{Y} should have been created by calling
  \code{\link{envelope}} with \code{savepatterns=TRUE}. The envelope of
  the summary function \code{fun} for these point patterns
  will be computed using the parameters specified in \code{\dots}.

  If \code{fun} is not given, then the new envelope will be computed
  from the same set of simulated function values
  that were used to compute \code{Y}.
  Also, \code{Y} should contain these function values,
  so \code{Y} should have been created by calling
  \code{\link{envelope}} with \code{savefuns=TRUE}. The envelope 
  will be computed using the parameters specified in \code{\dots}.
  
  See \code{\link{envelope}} for a full list of envelope parameters.
  Frequently-used parameters include \code{nrank} and \code{nsim} (to change the
  number of simulations used and the significance level of the
  envelope), \code{global} (to change from pointwise to global
  envelopes) and \code{VARIANCE} (to compute the envelopes from the sample
  moments instead of the ranks).
}
\value{
  An envelope (object of class \code{"envelope"}.
}
\seealso{
   \code{\link{envelope}}
}
\examples{
  data(cells)
  E <- envelope(cells, Kest, nsim=19, savefuns=TRUE, savepatterns=TRUE)
  E2 <- envelope(E, nrank=2)
  Eg <- envelope(E, global=TRUE)
  EG <- envelope(E, Gest)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
