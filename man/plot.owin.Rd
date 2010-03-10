\name{plot.owin}
\alias{plot.owin}
\title{Plot a Spatial Window}
\description{
  Plot a two-dimensional window of observation for a spatial point pattern
}
\usage{
 \method{plot}{owin}(x, main, add=FALSE, \dots, box, edge=0.04,
                      hatch=FALSE, angle=45, spacing=diameter(x)/50,
                      invert=FALSE)
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
  \item{hatch}{
    logical flag; if \code{TRUE}, the interior of the window will
    be shaded by a grid of parallel lines.
  }
  \item{angle}{
    orientation of the shading lines (in degrees anticlockwise
    from the \eqn{x} axis) when \code{hatch=TRUE}.
  }
  \item{spacing}{
    spacing between the shading lines,
    when \code{hatch=TRUE}.
  }
  \item{invert}{
    logical flag; when the window is a binary pixel mask,
    the mask colours will be inverted if \code{invert=TRUE}.
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

  When \code{x} is of type \code{"rectangle"} or \code{"polygonal"},
  it is plotted by the \R function \code{\link{polygon}}. To control the
  appearance (colour, fill density, line density etc) of the polygon
  plot, determine the required argument of \code{\link{polygon}}
  and pass it through \code{...} For example, to paint the
  interior of the polygon in red, use the argument
  \code{col="red"}. To draw the polygon edges in green,
  use \code{border="green"}. To suppress the drawing of polygon edges,
  use \code{border=NA}.

  When \code{x} is of type \code{"mask"}, it is plotted by
  \code{\link{image.default}}. The appearance of the image plot
  can be controlled by passing arguments to \code{\link{image.default}}
  through \code{...}. The default appearance can also be changed
  by setting the parameter \code{par.binary} of
  \code{\link{spatstat.options}}. 

  To zoom in (to view only a subset of the window at higher
  magnification), use the graphical arguments
  \code{xlim} and \code{ylim} to specify the desired rectangular field of
  view. (The actual field of view may be larger, depending on the
  graphics device).
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{plot.ppp}},
  \code{\link{polygon}},
  \code{\link{image.default}},
  \code{\link{spatstat.options}}
}
\examples{
  # rectangular window
   data(nztrees)
   plot(nztrees$window)
   abline(v=148, lty=2)

  # polygonal window
  data(demopat)
  w <- demopat$window
  plot(w)
  plot(w, col="red", border="green", lwd=2)
  plot(w, hatch=TRUE, lwd=2)

  # binary mask
  we <- erosion(w, 100, FALSE)
  plot(we)
  spatstat.options(par.binary=list(col=grey(c(0.5,1))))
  plot(we)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
