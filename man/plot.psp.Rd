\name{plot.psp}
\alias{plot.psp}
\title{plot a Spatial Line Segment Pattern}
\description{
  Plot a two-dimensional line segment pattern
}
\usage{
  \method{plot}{psp}(x, \dots, add=FALSE)
}
\arguments{
  \item{x}{
    The line segment pattern to be plotted.
    An object of class \code{"psp"},
    or data which can be converted into 
    this format by \code{\link{as.psp}()}.
  }
  \item{\dots}{
    extra arguments that will be passed to the plotting functions
    \code{\link{segments}} (to plot the segments)
    and 
    \code{\link{plot.owin}} (to plot the observation window).
  }
  \item{add}{
    Logical. If \code{TRUE}, the current plot is not erased;
    the segments are plotted on top of the
    current plot, and the window is not plotted.
  }
}
\value{
  \code{NULL}
}
\details{
  This is the \code{plot} method for 
  line segment pattern datasets (of class \code{"psp"},
  see \code{\link{psp.object}}).
  It plots both the observation window \code{x$window}
  and the line segments themselves.
  
  Plotting of the window \code{x$window} is performed by
  \code{\link{plot.owin}}. This plot may be modified 
  through the \code{...} arguments. 

  Plotting of the points themselves is performed
  by the standard R function \code{\link{segments}}.
  Its plotting behaviour may also be modified through the \code{...}
  arguments.

}
\seealso{
  \code{\link{psp.object}},
  \code{\link{plot}},
  \code{\link{par}},
  \code{\link{plot.owin}},
  \code{\link{symbols}}
}
\examples{
   a <- psp(runif(20), runif(20), runif(20), runif(20), window=owin())
   plot(a)
   plot(a, main="My title")
   plot(a, col="blue")
   plot(a, lwd=3)
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 \keyword{hplot}
