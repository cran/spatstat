\name{plot.owin}
\alias{plot.owin}
\title{Plot a Spatial Window}
\description{
  Plot a two-dimensional window of observation for a spatial point pattern
}
\usage{
 plot.owin(x, main, add=FALSE, \dots, box=TRUE, edge=0.04)
}
\arguments{
  \item{x}{
    The window to be plotted.
    An object of class \code{\link{owin}},
    or data which can be converted into 
    this format by \code{\link{as.owin}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{add}{
    logical flag: if \code{TRUE}, draw the window in 
    the current plot; if \code{FALSE}, generate a new plot.
  }
  \item{\dots}{
    extra arguments passed to the generic \code{\link{plot}} function.
  }
  \item{box}{
    logical flag; if \code{TRUE}, plot the enclosing rectangular box
  }
  \item{edge}{
    nonnegative number; the plotting region will have coordinate limits
    that are \code{1 + edge}
    times as large as the limits of the rectangular box
    that encloses the pattern.
  }
}
\value{
  none.
}
\details{
  This is the \code{plot} method for the class \code{\link{owin}}.
  The action is to plot the boundary of the window on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  If the window \code{x} is of type \code{"rectangle"} or \code{"polygonal"},
  the boundary of the window is plotted as a polygon or series of polygons.
  If \code{x} is of type \code{"mask"} 
  the discrete raster approximation of the window is displayed
  as a binary image (white inside the window, black outside).

  Graphical parameters controlling the display (e.g. setting the
  colours) may be passed directly via the \code{...} arguments,
  or indirectly reset using 
  \code{\link{spatstat.options}}.
  The parameter \code{par.binary} of \code{\link{spatstat.options}}
  controls the display when \code{x} is of type \code{"mask"}.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{plot.ppp}},
  \code{\link{spatstat.options}}
}
\examples{
\dontrun{
   library(spatstat)

  # rectangular window
   data(nztrees)
   plot(nztrees$window)
   abline(v=148, lty=2)

  # polygonal window
  plot(c(0,1),c(0,1),type="n")
  bdry <- locator()
  # click the vertices of a polygon (anticlockwise)
  w <- owin(c(0,1), c(0,1), poly=bdry)
  plot(w)

  # binary mask
  we <- erode.owin(w, 0.05, FALSE)
  plot(we)
  spatstat.options(par.binary=list(col=grey(c(0.5,1))))
  plot(we)
}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
