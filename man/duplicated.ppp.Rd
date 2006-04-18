\name{duplicated.ppp}
\alias{duplicated.ppp}
\title{Determine Duplicated Points in a Spatial Point Pattern}
\description{
  Determines which points in a spatial point pattern
  are duplicates of previous points, and returns a logical vector.
}
\usage{
 \method{duplicated}{ppp}(x, \dots)
}
\arguments{
  \item{x}{
    A spatial point pattern
    (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  A logical vector of length equal to the number of points in \code{x}.
}
\details{
  This is a method for the generic function \code{duplicated} for 
  point pattern datasets (of class \code{"ppp"}, see \code{\link{ppp.object}}).

  Two points in a point pattern are deemed to be identical
  if their \eqn{x,y} coordinates are the same,
  and their marks are also the same (if they carry marks).
  The Examples section illustrates how it is possible for
  a point pattern to contain a pair of identical points.

  This function determines which points in \code{x} duplicate
  other points that appeared earlier in the sequence. It 
  returns a logical vector with entries that are \code{TRUE}
  for duplicated points and \code{FALSE} for unique (non-duplicated)
  points.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{unique.ppp}},
  \code{\link{multiplicity.ppp}}
}
\examples{
   X <- ppp(c(1,1,0.5), c(2,2,1), window=square(3))
   duplicated(X)
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
