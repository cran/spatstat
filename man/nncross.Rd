\name{nncross}
\alias{nncross}
\title{Nearest Neighbour in Another Point Pattern}
\description{
  Given two point patterns \code{X} and \code{Y},
  finds the nearest neighbour in \code{Y} of each point of \code{X}.
}
\usage{
  nncross(X, Y)
}
\arguments{
  \item{X, Y}{Two point patterns (objects of class \code{"ppp"}).}
}
\details{
  Given two point patterns \code{X} and \code{Y} this
  function finds, for each point of \code{X}, the 
  the nearest point of \code{Y}. The distance between these points
  is also computed.

  The return value is a data frame, with rows corresponding to
  the points of \code{X}.  The first column gives the nearest neighbour
  distances (i.e. the \code{i}th entry is the distance 
  from the \code{i}th point of \code{X} to the nearest point of
  \code{Y}). The second column gives the indices of the nearest
  neighbours (i.e.\ the \code{i}th entry is the index of
  the nearest point in \code{Y}.)

  Note that this function is not symmetric in \code{X} and \code{Y}.
  To find the nearest neighbour in \code{X} of each point in \code{Y},
  just use \code{nncross(Y,X)}.
}
\value{
  A data frame with two columns:
  \item{dist}{Nearest neighbour distance}
  \item{which}{Nearest neighbour index in \code{Y}}
}
\examples{
  X <- runifpoint(15)
  Y <- runifpoint(20)
  N <- nncross(X,Y)$which
  # note that length(N) = 15
  plot(superimpose(X=X,Y=Y), main="nncross", cols=c("red","blue"))
  arrows(X$x, X$y, Y[N]$x, Y[N]$y, length=0.15)
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{math}
