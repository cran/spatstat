\name{pairwise.family}
\alias{pairwise.family}
\title{Pairwise Interaction Process Family}
\description{
  An object describing the family of all pairwise interaction point processes
}
\details{
  \bold{Advanced Use Only!}
  
  This structure would not normally be touched by
  the user. It describes 
  the pairwise interaction family of point process models.
 
  If you need to create a specific pairwise interaction model for use in 
  modelling, use the function \code{\link{Pairwise}} or one of the existing
  functions (listed below, starting with a capital letter).
 
  Anyway, \code{pairwise.family} is an object of class \code{"isf"}
  containing a function \code{pairwise$eval} for
  evaluating the sufficient statistics of any pairwise interaction
  point process model taking an exponential family form. 
} 
\seealso{
  \code{\link{pairsat.family}},
  \code{\link{ord.family}},
  \code{\link{Poisson}},
  \code{\link{Pairwise}},
  \code{\link{PairPiece}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{Geyer}},
  \code{\link{Saturated}},
  \code{\link{Ord}},
  \code{\link{OrdThresh}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
