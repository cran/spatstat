\name{gridweights}
\alias{gridweights}
\title{Compute Quadrature Weights Based on Grid Counts}
\description{
  Computes quadrature weights for a given set of points,
  using the ``counting weights'' for a grid of rectangular tiles.
}
\usage{
 gridweights(X, nx, ny, window=NULL)
}
\arguments{
  \item{X}{Data defining a point pattern.}
  \item{nx}{Number of tiles
    in each row of the rectangular grid.
  }
  \item{ny}{Number of tiles
    in each column of the rectangular grid.
  }
  \item{window}{Default window for the point pattern}
}
\value{
  Vector of nonnegative weights for each point in \code{X}.
}
\details{
  This function computes a set of quadrature weights
  for a given pattern of points
  (typically comprising both ``data'' and `dummy'' points).
  See \code{\link{quad.object}} for an explanation of quadrature
  weights and quadrature schemes.

  The weights are computed by the ``counting weights'' rule
  based on a regular grid of rectangular tiles.
  First \code{X} and (optionally) \code{window} are converted into a
  point pattern object. Then the bounding rectangle of the window of
  the point pattern is
  divided into a regular \eqn{nx \times ny}{nx * ny} grid of rectangular tiles.
  The weight attached to a point of \code{X} is the area of the tile
  in which it lies, divided by the number of points of \code{X} lying in
  that tile.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{dirichlet.weights}}
}
\examples{
  \dontrun{
  library(spatstat)
  }
  Q <- quadscheme(runifpoispp(10))
  X <- as.ppp(Q) # data and dummy points together
  w <- gridweights(X, 10, 10)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
