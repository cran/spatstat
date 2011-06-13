\name{methods.distfun}
\alias{methods.distfun}
\alias{contour.distfun}
\alias{persp.distfun}
\alias{plot.distfun}
\alias{print.distfun}
\Rdversion{1.1}
\title{
  Methods for Distance Functions
}
\description{
  Methods for objects of the class \code{"distfun"}.
}
\usage{
\method{contour}{distfun}(x, \dots)
\method{persp}{distfun}(x, \dots)
\method{plot}{distfun}(x, \dots)
\method{print}{distfun}(x, \dots)
}
\arguments{
  \item{x}{
    Object of class \code{"distfun"} representing a
    distance function.
  }
  \item{\dots}{
    Named arguments controlling the plot. See Details.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}},
  \code{\link{plot}},
  \code{\link{contour}}
  and \code{\link{persp}}
  for the class \code{"distfun"} of distance functions.
  See \code{\link{distfun}} for explanation about this class.
  
  The \code{print} method prints a description of the domain of the
  distance function and the spatial object of which it is the distance
  function. Any additional arguments \code{\dots} are ignored.
  
  The \code{plot}, \code{contour} and \code{persp} methods first convert
  \code{x} to a pixel image object using \code{\link{as.im}},
  then display it using \code{\link{plot.im}}, \code{\link{contour.im}} or
  \code{\link{persp.im}}.
  Additional arguments \code{\dots} are either
  passed to \code{\link{as.im.function}} to
  control the spatial resolution of the pixel image, or passed to 
  \code{\link{contour.im}},
  \code{\link{persp.im}} or 
  \code{\link{plot.im}} to control the appearance of the plot. 
}
\value{
  \code{NULL}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{distfun}},
  \code{\link{as.im}},
  \code{\link{plot.im}},
  \code{\link{persp.im}},
  \code{\link{contour.im}},
  \code{\link{spatstat.options}}
}
\examples{
   data(letterR)
   f <- distfun(letterR)
   contour(f)
   contour(f, W=owin(c(1,5),c(-1,4)), eps=0.1)
}
\keyword{spatial}
\keyword{methods}
