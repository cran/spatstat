\name{rThomas}
\alias{rThomas}
\title{Simulate Thomas Process}
\description{
  Generate a random point pattern using the
  Thomas cluster process.
}
\usage{
 rThomas(lambda, sigma, mu, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process of cluster centres.
    A single positive number.
  }
  \item{sigma}{
    Standard deviation of displacement of a point from its cluster centre.
  }
  \item{mu}{
    Expected number of points per cluster.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of the
  Thomas process, a special case of the Neyman-Scott process.
  
  The algorithm 
  generates a uniform Poisson point process of ``parent'' points 
  with intensity \code{lambda}. Then each parent point is
  replaced by a random cluster of points, the number of points
  per cluster being Poisson (\code{mu}) distributed, and their
  positions being isotropic Gaussian displacements from the
  cluster parent location.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rNeymanScott}}
}
\examples{
  X <- rThomas(10, 0.2, 5)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
