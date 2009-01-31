\name{runifpointOnLines}
\alias{runifpointOnLines}
\title{Generate N Uniform Random Points On Line Segments}
\description{
  Given a line segment pattern, generate a random point pattern
  consisting of \code{n} points uniformly distributed on the
  line segments.
}
\usage{
runifpointOnLines(n, L)
}
\arguments{
  \item{n}{Number of points to generate.}
  \item{L}{Line segment pattern (object of class \code{"psp"})
    on which the points should lie.
  }
}
\details{
  This command generates a point pattern consisting of \code{n}
  independent random points, each point uniformly distributed on
  the line segment pattern. This means that, for each random point,
  \itemize{
    \item the
    probability of falling on a particular segment is proportional to the
    length of the segment; and
    \item given that the point falls on a particular segment,
    it has uniform probability density along that segment.
  }
}
\value{
  A point pattern (object of class \code{"ppp"}) with the same
  window as \code{L}.
}
\seealso{
  \code{\link{psp}},
  \code{\link{ppp}},
  \code{\link{pointsOnLines}},
  \code{\link{runifpoint}}
}
\examples{
  X <- psp(runif(10), runif(10), runif(10), runif(10),  window=owin())
  Y <- runifpointOnLines(20, X)
  plot(X, main="")
  plot(Y, add=TRUE)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
