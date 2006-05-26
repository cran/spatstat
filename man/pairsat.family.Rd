\name{pairsat.family}
\alias{pairsat.family}
\title{Saturated Pairwise Interaction Point Process Family}
\description{
  An object describing the Saturated Pairwise Interaction
  family of point process models
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes the ``saturated pairwise interaction''
  family of point process models.
 
  If you need to create a specific interaction model for use in 
  spatial pattern analysis, use the function \code{\link{Saturated}()}
  or the one existing implementation of a model in this family,
  \code{\link{Geyer}()}.
 
  Geyer (1999) introduced the ``saturation process'', a modification of the
  Strauss process in which the total contribution
  to the potential  from each point (from its pairwise interaction with all 
  other points) is trimmed to a maximum value \eqn{c}. 
  This model is implemented in the function \code{\link{Geyer}()}.
 
  The present class \code{pairsat.family} is the 
  extension of this saturation idea to all pairwise interactions.
  Note that the resulting models are no longer pairwise interaction
  processes - they have interactions of infinite order.
 
  \code{pairsat.family} is an object of class \code{"isf"} 
  containing a function \code{pairwise$eval} for
  evaluating the sufficient statistics of any saturated pairwise interaction
  point process model in which the original pair potentials 
  take an exponential family form. 
} 
\references{
Geyer, C.J. (1999)
Likelihood Inference for Spatial Point Processes.
Chapter 3 in 
O.E. Barndorff-Nielsen, W.S. Kendall and M.N.M. Van Lieshout (eds)
\emph{Stochastic Geometry: Likelihood and Computation},
Chapman and Hall / CRC, 
Monographs on Statistics and Applied Probability, number 80.
Pages 79--140.
}
\seealso{
  \code{\link{pairwise.family}},
  \code{\link{Geyer}}
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
