\name{clickpoly}
\alias{clickpoly}
\title{Interactively Define a Polygon}
\description{
  Allows the user to create a polygon by
  point-and-click in the display.
}
\usage{
  clickpoly(add=FALSE, n=NULL, ...)
}
\arguments{
  \item{add}{
    Logical value indicating whether to create a new plot
    (\code{add=FALSE}) or draw over the existing plot (\code{add=TRUE}).
  }
  \item{n}{
    Number of vertices of the polygon (if this is predetermined).
  }
  \item{\dots}{
    Arguments passed to \code{\link{locator}} to control the
    interactive plot.
  }
}
\value{
  A window (object of class \code{"owin"}) representing the polygon.
}
\details{
  This function allows the user to create a polygonal window 
  by interactively clicking on the screen display.

  The user is prompted to point the mouse at any desired locations
  for the polygon vertices,
  and click the left mouse button to add each point.
  Interactive input stops after \code{n} clicks (if \code{n} was given)
  or when the middle mouse button is pressed.

  The return value is a window (object of class \code{"owin"})
  representing the polygon.

  This function uses the \R{} command \code{\link{locator}} to
  input the mouse clicks. It only works on screen devices such as
  \sQuote{X11}, \sQuote{windows} and \sQuote{quartz}. Arguments that can be
  passed to \code{\link{locator}} through \code{\dots} include
  \code{pch} (plotting character), \code{cex} (character expansion
  factor) and \code{col} (colour). See \code{\link{locator}}
  and \code{\link{par}}.
}
\seealso{
  \code{\link{identify.ppp}},
  \code{\link{clickppp}},
  \code{\link{locator}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{iplot}
