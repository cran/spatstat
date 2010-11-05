\name{rpoisppOnLines}
\alias{rpoisppOnLines}
\title{Generate Poisson Point Pattern on Line Segments}
\description{
  Given a line segment pattern, generate a Poisson random point pattern
  on the line segments.
}
\usage{
rpoisppOnLines(lambda, L, lmax = NULL, ...)
}
\arguments{
  \item{lambda}{Intensity of the Poisson process.
    A single number, a \code{function(x,y)}, or a pixel image
    (object of class \code{"im"}).
  }
  \item{L}{Line segment pattern (object of class \code{"psp"})
    on which the points should be generated.
  }
  \item{lmax}{Maximum possible value of \code{lambda}
    if it is a function or a pixel image.
  }
  \item{\dots}{Additional arguments passed to \code{lambda} if it is a
    function.
  }
}
\details{
  This command generates a Poisson point process on the one-dimensional
  system of line segments in \code{L}. The result is a point pattern
  consisting of points lying on the line segments in \code{L}.
  The number of random points falling on any given line segment follows
  a Poisson distribution. The patterns of points on different segments
  are independent.

  The intensity \code{lambda} is the expected number of points
  per unit \bold{length} of line segment. 
  It may be constant, or it may depend on spatial location.

  The argument \code{lambda} may be a single number,
  or a \code{function(x,y)}, or a pixel image (object of class
  \code{"im"}). In the two latter cases, the rejection method is used.

  The rejection method requires knowledge of \code{lmax},
  the maximum possible value of \code{lambda}. If \code{lmax} is not
  given, it will be computed approximately, by sampling many values of
  \code{lambda}. 
}
\value{
  A point pattern (object of class \code{"ppp"}) in the same window
  as \code{L}.
}
\seealso{
  \code{\link{psp}},
  \code{\link{ppp}},
  \code{\link{runifpointOnLines}},
  \code{\link{rpoispp}}
}
\examples{
  L <- psp(runif(10), runif(10), runif(10), runif(10),  window=owin())
  plot(L, main="")

  # uniform intensity
  Y <- rpoisppOnLines(4, L)
  plot(Y, add=TRUE, pch="+")

  # intensity is a function
  Y <- rpoisppOnLines(function(x,y){ 10 * x^2}, L, 10)
  plot(L, main="")
  plot(Y, add=TRUE, pch="+")

  # intensity is an image
  Z <- as.im(function(x,y){10 * sqrt(x+y)}, unit.square())
  Y <- rpoisppOnLines(Z, L, 10)
  plot(L, main="")
  plot(Y, add=TRUE, pch="+")
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
