\name{plot.im}
\alias{plot.im}
\alias{image.im}
\title{Plot a Pixel Image}
\description{
  Plot a pixel image.
}
\usage{
 plot.im(x, \dots, ribbon=TRUE, ribsep=0.15, ribwid=0.05, ribn=1024)
}
\arguments{
  \item{x}{
    The pixel image to be plotted.
    An object of class \code{"im"} (see \code{\link{im.object}}).
  }
  \item{\dots}{
    Extra arguments passed to \code{\link{image}} to control the plot.
  }
  \item{ribbon}{
    Logical flag indicating whether to display a ribbon
    showing the colour map.
  }
  \item{ribsep}{
    Factor controlling the space between the ribbon and the image.
  }
  \item{ribwid}{
    Factor controlling the width of the ribbon.
  }
  \item{ribn}{
    Number of different values to display in the ribbon.
  }
}
\value{
  none.
}
\details{
  This is the \code{plot} method for the class \code{"im"}.
  [It is also the \code{image} method for \code{"im"}.]

  The pixel image \code{x} is displayed on the current plot device,
  using equal scales on the \code{x} and \code{y} axes.

  If \code{ribbon=TRUE}, a legend will be plotted at the right of the
  image. The legend consists of a colour ribbon and an axis with tick-marks,
  showing the correspondence between the pixel values and the colour map.

  Arguments \code{ribsep, ribwid, ribn} control the appearance of the
  ribbon.
  The width of the ribbon is \code{ribwid} times the size of the pixel
  image, where `size' means the larger of the width and the height.
  The distance separating the ribbon and the image is \code{ribsep} times
  the size of the pixel image. The ribbon contains \code{ribn}
  different numerical values, evenly spaced between the minimum and
  maximum pixel values in the image \code{x}, rendered according to
  the chosen colour map.
  
  Graphical parameters controlling the display of both the pixel image
  and the ribbon are passed through the \code{...} arguments
  directly to the function \code{\link{image.default}}.
  In particular, the colour map is controlled by the argument
  \code{col}. See the examples in \code{\link{image.default}}.
}
\seealso{
  \code{\link{im.object}}
}
\examples{
   # an image
   Z <- setcov(owin())
   plot(Z)
   plot(Z, col=terrain.colors(128), axes=FALSE)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
