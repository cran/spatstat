\name{rNeymanScott}
\alias{rNeymanScott}
\title{Simulate Neyman-Scott Process}
\description{
  Generate a random point pattern using the
  Neyman-Scott cluster process.
}
\usage{
 rNeymanScott(lambda, rmax, rcluster, win = owin(c(0,1),c(0,1)), \dots)
}
\arguments{
  \item{lambda}{
    Intensity of the Poisson process of cluster centres.
    A single positive number.
  }
  \item{rmax}{
    Maximum radius of a random cluster.
  }
  \item{rcluster}{
    A function which generates random clusters.
  }
  \item{win}{
    Window in which to simulate the pattern.
    An object of class \code{"owin"}
    or something acceptable to \code{\link{as.owin}}.
  }
  \item{\dots}{
    Arguments passed to \code{rcluster}
  }
}
\value{
  The simulated point pattern (an object of class \code{"ppp"}).
}
\details{
  This algorithm generates a realisation of the
  general Neyman-Scott process, with the cluster mechanism
  given by the function \code{rcluster}. 
  The clusters must have a finite maximum possible radius \code{rmax}.

  We algorithm 
  generates a uniform Poisson point process of ``parent'' points 
  with intensity \code{lambda}. Then each parent point is
  replaced by a random cluster of points, created by calling the
  function \code{rcluster}.

  The function \code{rcluster} should expect to be called as
  \code{rcluster(xp[i],yp[i],\dots)} for each parent point at a location
  \code{(xp[i],yp[i])}. The return value of \code{rcluster}
  should be a list with elements
  \code{x,y} which are vectors of equal length giving the absolute
  \eqn{x} and \code{y} coordinates of the points in the cluster.
}
\seealso{
\code{\link{rpoispp}},
\code{\link{rMatClust}}
}
\examples{
  \dontrun{library(spatstat)}
  nclust <-  function(x0, y0, radius, n) {
                           return(runifdisc(n, radius, x0, y0))
                         }
  X <- rNeymanScott(10, 0.2, nclust, radius=0.2, n=5)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
