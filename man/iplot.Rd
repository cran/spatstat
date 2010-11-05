\name{iplot}
\alias{iplot}
\title{Point and Click Interface for Displaying a Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
  with interactive (point-and-click) control over the plot.
}
\usage{
 iplot(x, xname)
}
\arguments{
  \item{x}{
    The spatial point pattern to be plotted.
    An object of class \code{"ppp"}.
  }
  \item{xname}{
    Optional. Character string to use as the title of the dataset.
  }
}
\value{
  \code{NULL}.
}
\details{
  This function generates a plot of a spatial point pattern dataset
  (object of class \code{"ppp"}) and allows 
  interactive control over the appearance of the plot
  using a point-and-click interface.

  A new popup window is launched.
  The point pattern \code{x} is displayed in the left half of the window
  using \code{\link{plot.ppp}} with the default values of all
  the plot parameters. The right side of the window contains
  buttons and sliders allowing the user to change the header text,
  the plot symbols, the scale used to represent numeric marks,
  and so on. For a multitype point pattern, the user
  can also switch between viewing all points in a single display,
  and splitting the points into separate patterns according to type.
}
\seealso{
  \code{\link{istat}},
  \code{\link{ppp.object}},
  \code{\link{plot.ppp}}
}
\examples{
   if(interactive()) {
      data(cells)
      iplot(cells)
   }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{hplot}
