\name{runifpoint}
\alias{runifpoint}
\title{Generate N Uniform Random Points}
\description{
  Generate a random point pattern
  containing \eqn{n} independent uniform random points.
}
\usage{
 runifpoint(n, win=owin(c(0,1),c(0,1)), giveup=1000)
}
\arguments{
  \item{n}{
    Number of points.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{giveup}{
    Number of attempts in the rejection method after which
    the algorithm should stop trying to generate new points.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  If \code{win} is a rectangle then this algorithm simply
  generates \eqn{n} independent random points uniformly distributed
  in the rectangle, by assigning uniform random values to their
  cartesian coordinates.

  If \code{win} is not a rectangle, the algorithm uses the rejection
  method. It finds a rectangle enclosing the window,
  generates points in this rectangle, and tests whether they fall in
  the desired window. It gives up after \code{giveup} consecutive
  failures have occurred.
}
\seealso{
\code{\link{ppp.object}},
\code{\link{owin.object}}
}
\examples{
 # 100 random points in the unit square
 pp <- runifpoint(100)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
