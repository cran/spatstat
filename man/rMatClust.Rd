\name{rMatClust}
\alias{rMatClust}
\title{Simulate Matern Cluster Process}
\description{
  Generate a random point pattern, a simulated realisation of the
  Mat\'ern Cluster Process.
}
\usage{
 rMatClust(kappa, r, mu, win = owin(c(0,1),c(0,1)))
}
\arguments{
  \item{kappa}{
    Intensity of the Poisson process of cluster centres.
    A single positive number, a function, or a pixel image.
  }
  \item{r}{
    Radius parameter of the clusters.
  }
  \item{mu}{
    Mean number of points per cluster.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).

  Additionally,  some intermediate results of the simulation are
  returned as attributes of this point pattern.
  See \code{\link{rNeymanScott}}.
}
\details{
  This algorithm generates a realisation of Mat\'ern's cluster process
  inside the window \code{win}. The process is constructed by first
  generating a Poisson point process of ``parent'' points 
  with intensity \code{kappa}. Then each parent point is
  replaced by a random cluster of points, the number of points in each
  cluster being random with a Poisson (\code{mu}) distribution,
  and the points being placed independently and uniformly inside
  a disc of radius \code{r} centred on the parent point.

  In this implementation, parent points are not restricted to lie in the
  window; the parent process is effectively the uniform
  Poisson process on the infinite plane.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rNeymanScott}}
}
\examples{
 pp <- rMatClust(10, 0.05, 4)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{datagen}

