\name{applynbd}
\alias{applynbd}
\title{Apply Function to Every Neighbourhood in a Point Pattern}
\description{
  Visit each point in a point pattern, find the neighbouring points,
  and apply a given function to them.
}
\usage{
   applynbd(X, FUN, N, R, criterion, exclude=FALSE, \dots)
}
\arguments{
  \item{X}{
    Point pattern. 
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}}.
  }
  \item{FUN}{
    Function to be applied to each neighbourhood.
    The arguments of \code{FUN} are described under \bold{Details}.
  }
  \item{N}{
    Integer. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of the
    \code{N} points of \code{X} which are closest to it.
    This argument is incompatible with \code{R} and \code{criterion}.
  }
  \item{R}{
    Nonnegative numeric value. If this argument is present,
    the neighbourhood of a point of \code{X} is defined to consist of
    all points of \code{X} which lie within a distance \code{R}
    of it.
    This argument is incompatible with \code{N} and \code{criterion}.
  }
  \item{criterion}{
    Function. If this argument is present, 
    the neighbourhood of a point of \code{X} is determined by
    evaluating this function. See under \bold{Details}.
    This argument is incompatible with \code{N} and \code{R}.
  }
  \item{exclude}{
    Logical. If \code{TRUE} then the point currently being visited
    is excluded from its own neighbourhood.
  }
  \item{\dots}{
    extra arguments passed to the function \code{FUN}.
    They must be given in the form \code{name=value}.
  }
}
\value{
  Similar to the result of \code{\link{apply}}.
  If each call to \code{FUN} returns a single numeric value,
  the result is a vector of dimension \code{X$n}, the number of points
  in \code{X}.
  If each call to \code{FUN} returns a vector of the same length
  \code{m}, then the result is a matrix of dimensions \code{c(m,n)};
  note the transposition of the indices, as usual for the family of
  \code{apply} functions.
  If the calls to \code{FUN} return vectors of different lengths,
  the result is a list of length \code{X$n}. 
}
\details{
  This is an analogue of \code{\link{apply}}
  for point patterns. It visits each point in the point pattern \code{X},
  determines which points of \code{X} are ``neighbours'' of the current
  point, applies the function \code{FUN} to this neighbourhood,
  and collects the values returned by \code{FUN}.

  The definition of ``neighbours'' depends on the arguments
  \code{N}, \code{R} and \code{criterion}, exactly one of which
  must be given.
  Also the argument \code{exclude} determines whether
  the current point is excluded from its own neighbourhood.

  If \code{N} is given, then the neighbours of the current
  point are the \code{N} points of \code{X} which are closest to
  the current point (including the current point itself
  unless \code{exclude=TRUE}).
  If \code{R} is given, then the neighbourhood of the current point
  consists of all points of \code{X} which lie closer than a distance \code{R}
  from the current point.
  If \code{criterion} is given, then it must be a function
  with two arguments \code{dist} and \code{drank} which will be
  vectors of equal length.
  The interpretation is that \code{dist[i]} will be the
  distance of a point from the current point, and
  \code{drank[i]} will be the rank of that distance (the three points
  closest to the current point will have rank 1, 2 and 3).
  This function must return a logical vector of the same length
  as \code{dist} and \code{drank} whose \code{i}-th entry is
  \code{TRUE} if the corresponding point should be included in
  the neighbourhood. See the examples below.

  Each point of \code{X} is visited; the neighbourhood
  of the current point is determined, and stored as a point pattern
  \code{Y}; then the function \code{FUN} is called as
  \code{FUN(Y, current, dists, dranks, \dots)}
  where \code{current} is a \code{list(x,y)} giving the location of the
  current point, \code{dists} is a vector of distances from the current
  point to each of the points in \code{Y}, 
  \code{dranks} is a vector of the ranks of these distances
  with respect to the full point pattern \code{X},
  and \code{\dots} are the arguments passed from the call to
  \code{applynbd}.

  The results of each call to \code{FUN} are collected and returned
  according to the usual rules for \code{\link{apply}} and its
  relatives. See \bold{Value} above.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{apply}}
}
\examples{
  data(redwood)

  # count the number of points within radius 0.2 of each point of X
  nneighbours <- applynbd(redwood, R=0.2, function(Y, ...){Y$n-1})
  # equivalent to:
  nneighbours <- applynbd(redwood, R=0.2, function(Y, ...){Y$n}, exclude=TRUE)

  # compute the distance to the second nearest neighbour of each point
  secondnndist <- applynbd(redwood, N = 2, function(Y, cur, d, r){max(d)}, exclude=TRUE)

  # marked point pattern
  data(longleaf)
  \testonly{
	# smaller dataset
	longleaf <- longleaf[seq(1, longleaf$n, by=80)]
  }
  # compute the median of the marks of all neighbours of a point
  dbh.med <- applynbd(longleaf, R=40, exclude=TRUE,
                 function(Y,...) { median(Y$marks)})


  # ANIMATION explaining the definition of the K function
  # (arguments `fullpicture' and 'rad' are passed to FUN)

  showoffK <- function(Y, u, d, r, fullpicture,rad) { 
	plot(fullpicture, main="")
	points(Y[r>0], cex=2)
	points(u[1],u[2],pch="+",cex=3)
	theta <- seq(0,2*pi,length=100)
	polygon(u[1]+ rad * cos(theta),u[2]+rad*sin(theta))
	text(u[1]+rad/3,u[2]+rad/2,Y$n-1,cex=3)
	if(runif(1) < 0.1)
		system("sleep 1")
	return(Y$n - 1)
  }
  applynbd(redwood, R=0.2, showoffK, fullpicture=redwood, rad=0.2)

  # animation explaining the definition of the G function

  showoffG <- function(Y, u, d, r, fullpicture) { 
	plot(fullpicture, main="")
	points(Y, cex=2)
	points(u[1],u[2],pch="+",cex=3)
	v <- c(Y$x[1],Y$y[1])
	segments(u[1],u[2],v[1],v[2],lwd=2)
	w <- (u + v)/2
	nnd <- sqrt(sum((u-v)^2))
	text(w[1],w[2],round(nnd,3),cex=2)
	if(runif(1) < 0.1)
		system("sleep 1")
	return(nnd)
  }

  data(cells)
  applynbd(cells, N=1, showoffG, exclude=TRUE, fullpicture=cells)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
