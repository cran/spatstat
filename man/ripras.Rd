\name{ripras}
\alias{ripras}
\title{Estimate window from points alone}
\description{
  Given an observed pattern of points,
  computes the Ripley-Rasson estimate of 
  the spatial domain from which they came.
}
\usage{
 ripras(x, y=NULL)
}
\arguments{
  \item{x}{vector of \code{x} coordinates of observed points,
	   or a 2-column matrix giving \code{x,y} coordinates,
	   or a list with components \code{x,y} giving coordinates.}
  \item{y}{(optional) vector of \code{y} coordinates of observed points,
           if \code{x} is a vector.}
}
\value{
  A window (an object of class \code{"owin"}).
}
\details{
  Given an observed pattern of points with coordinates 
  given by \code{x} and \code{y}, this function computes 
  an estimate due to Ripley and Rasson (1977) of the
  spatial domain from which the points came. 

  The points are
  assumed to have been generated independently and uniformly
  distributed inside an unknown domain \eqn{D}. The maximum
  likelihood estimate of \eqn{D} is the convex hull of the 
  points. Analogously to the problems of estimating the endpoint
  of a uniform distribution, the MLE is not optimal.
  Ripley and Rasson's estimator is a rescaled copy of the convex hull,
  centred at the centroid of the convex hull.
  The scaling factor is 
  \eqn{1/sqrt(1 - m/n)}{1/\sqrt{1 - \frac m n}}
  where \eqn{n} is the number of data points and 
  \eqn{m} the number of vertices of the convex hull.
}
\seealso{
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  \dontrun{library(spatstat)
  plot(owin())
  }
  x <- runif(30)
  y <- runif(30)
  \dontrun{points(x,y)}
  w <- ripras(x,y)
  \dontrun{
  plot(w, box=FALSE)
  points(x,y)
  }
}
\references{
  Ripley, B.D. and Rasson, J.-P. (1977)
  Finding the edge of a Poisson forest.
  \emph{Journal of Applied Probability},
  \bold{14}, 483 -- 491.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}




