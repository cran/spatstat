\name{nncross}
\alias{nncross}
\title{Nearest Neighbours Between Two Patterns}
\description{
  Given two point patterns \code{X} and \code{Y},
  finds the nearest neighbour in \code{Y} of each point of \code{X}.
  Alternatively \code{Y} may be a line segment pattern.
}
\usage{
  nncross(X, Y, iX=NULL, iY=NULL)
}
\arguments{
  \item{X}{Point pattern (object of class \code{"ppp"}).}
  \item{Y}{Either a point pattern (object of class \code{"ppp"})
    or a line segment pattern (object of class \code{"psp"}).}
  \item{iX, iY}{Optional identifiers, applicable only in the case where
    \code{Y} is a point pattern, used to determine whether a point in
    \code{X} is identical to a point in \code{Y}. See Details}
}
\details{
  Given two point patterns \code{X} and \code{Y} this
  function finds, for each point of \code{X}, 
  the nearest point of \code{Y}. The distance between these points
  is also computed.

  Alternatively if \code{X} is a point pattern and \code{Y} is a line
  segment pattern, the function finds the nearest line segment to each point
  of \code{X}, and computes the distance.

  The return value is a data frame, with rows corresponding to
  the points of \code{X}.  The first column gives the nearest neighbour
  distances (i.e. the \code{i}th entry is the distance 
  from the \code{i}th point of \code{X} to the nearest element of
  \code{Y}). The second column gives the indices of the nearest
  neighbours (i.e.\ the \code{i}th entry is the index of
  the nearest element in \code{Y}.)

  Note that this function is not symmetric in \code{X} and \code{Y}.
  To find the nearest neighbour in \code{X} of each point in \code{Y},
  where \code{Y} is a point pattern, use \code{nncross(Y,X)}.

  The arguments \code{iX} and \code{iY} are used when
  the two point patterns \code{X} and \code{Y} have some points in
  common.  In this situation \code{nncross(X, Y)} would return some zero
  distances. To avoid this, attach a unique integer identifier to
  each point, such that two points are identical if their
  identifying numbers are equal. Let \code{iX} be the vector of
  identifier values for the points in \code{X}, and \code{iY}
  the vector of identifiers for points in \code{Y}. Then the code
  will only compare two points if they have different values of the
  identifier. See the Examples.
}
\value{
  A data frame with two columns:
  \item{dist}{Nearest neighbour distance}
  \item{which}{Nearest neighbour index in \code{Y}}
}
\seealso{
  \code{\link{nndist}} for nearest neighbour
  distances in a single point pattern.
}
\examples{
  # two different point patterns
  X <- runifpoint(15)
  Y <- runifpoint(20)
  N <- nncross(X,Y)$which
  # note that length(N) = 15
  plot(superimpose(X=X,Y=Y), main="nncross", cols=c("red","blue"))
  arrows(X$x, X$y, Y[N]$x, Y[N]$y, length=0.15)

  # two patterns with some points in common
  Z <- runifpoint(50)
  X <- Z[1:30]
  Y <- Z[20:50]
  iX <- 1:30
  iY <- 20:50
  N <- nncross(X,Y, iX, iY)$which
  plot(superimpose(X=X, Y=Y), main="nncross", cols=c("red","blue"))
  arrows(X$x, X$y, Y[N]$x, Y[N]$y, length=0.15)

  # point pattern and line segment pattern
  X <- runifpoint(15)
  Y <- rpoisline(10)
  N <- nncross(X,Y)
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
