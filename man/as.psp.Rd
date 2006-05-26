\name{as.psp}
\alias{as.psp}
\alias{as.psp.psp}
\alias{as.psp.data.frame}
\alias{as.psp.matrix}
\alias{as.psp.owin}
\alias{as.psp.default}
\title{Convert Data To Class psp}
\description{
  Tries to coerce any reasonable kind of data to a line segment pattern
  (an object of class \code{"psp"})
  for use by the \pkg{spatstat} package.
}
\synopsis{
  as.psp(x, \dots)
  as.psp.psp(x, \dots, fatal=TRUE)
  as.psp.data.frame(x, \dots, window=NULL, marks=NULL, fatal=TRUE)
  as.psp.matrix(x, \dots, window=NULL, marks=NULL, fatal=TRUE)
  as.psp.owin(x, \dots, fatal=TRUE)
  as.psp.default(x, \dots, window=NULL, marks=NULL, fatal=TRUE)
}
\usage{
  as.psp(x, \dots, window=, marks=, fatal=TRUE)
}
\arguments{
  \item{x}{Data which will be converted into a line segment pattern}
  \item{window}{Data which define a window for the pattern
    when \code{x} does not contain a window}
  \item{\dots}{Ignored.}
  \item{marks}{(Optional) vector of marks for the pattern}
  \item{fatal}{Logical value. See Details.}
}
\value{
  An object of class \code{"psp"} (see \code{\link{psp.object}})
  describing the line segment pattern and its window of observation.
  The value \code{NULL} may also be returned; see Details.
}
\details{
  Converts the dataset \code{x} to a line segment pattern
  (an object of class \code{"psp"}; see \code{\link{psp.object}} for
  an overview).

  This function is normally used to convert an existing line segment pattern
  dataset, stored in another format, to the \code{"psp"} format.
  To create a new point pattern from raw data such as \eqn{x,y}
  coordinates, it is normally easier to use the creator function
  \code{\link{psp}}.

  The dataset \code{x} may be:
  \itemize{
    \item
    an object of class \code{"psp"}  
    \item
    a data frame with columns named \code{x0, y0, x1, y1}
    that will be interpreted as the coordinates of the endpoints
    of the segments
    \item 
    a data frame or matrix with exactly 4 columns
    that will be interpreted as the coordinates \code{x0, y0, x1, y1}
    of the endpoints of the segments
    \item
    a data frame with columns named \code{xmid, ymid, length, angle}
    that will be interpreted as the coordinates of the segment
    midpoints, the lengths of the segments, and the orientations of the
    segments in radians
    \item
    a structure with elements named  \code{x0, y0, x1, y1}
    or elements named \code{xmid, ymid, length, angle}
    which will be interpreted as above.
    \item
    an object of class \code{"owin"} representing a spatial window;
    it must be of type \code{"rectangle"} or \code{"polygonal"}.
    The boundary edges of the window will be extracted as a line segment
    pattern.
  }
  The argument \code{window} is converted to a window object
  by the function \code{\link{as.owin}}.

  The argument \code{fatal} indicates what to do when
  the data cannot be converted to a line segment pattern.
  If \code{fatal=TRUE}, a fatal error
  will be generated; if \code{fatal=FALSE}, the
  value \code{NULL} is returned.

  The function \code{as.psp} is generic, with methods for the
  classes \code{"psp"}, \code{"data.frame"}, \code{"matrix"}
  and a default method.
  
  Point pattern datasets can also be created by the function
  \code{\link{psp}}.
}
\seealso{
  \code{\link{psp}},
  \code{\link{psp.object}},
  \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\examples{
   mat <- matrix(runif(40), ncol=4)
   a <- as.psp(mat, window=owin())
   stuff <- list(xmid=runif(10),
                 ymid=runif(10),
                 length=rep(0.1, 10),
                 angle=runif(10, 0, 2 * pi))
   a <- as.psp(stuff, window=owin())
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{manip}
