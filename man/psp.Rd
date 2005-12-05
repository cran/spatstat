\name{psp}
\alias{psp}
\title{Create a Line Segment Pattern}
\description{
  Creates an object of class \code{"psp"} representing 
  a line segment pattern in the two-dimensional plane.
}
\usage{
  psp(x0,y0, x1, y1, window, marks=NULL)
}
\arguments{
  \item{x0}{Vector of \eqn{x} coordinates of first endpoint of each segment}
  \item{y0}{Vector of \eqn{y} coordinates of first endpoint of each segment}
  \item{x1}{Vector of \eqn{x} coordinates of second endpoint of each segment}
  \item{y1}{Vector of \eqn{y} coordinates of second endpoint of each segment}
  \item{window}{window of observation,
    an object of class \code{"owin"}}
  \item{marks}{(optional) vector of mark values}
}
\value{
  An object of class \code{"psp"} 
  describing a line segment pattern in the two-dimensional plane
  (see \code{\link{psp.object}}).
}
\details{
  In the \code{spatstat} library, a spatial pattern of line segments is
  described by an object of class \code{"psp"}. This function
  creates such objects.

  The vectors \code{x0}, \code{y0}, \code{x1} and \code{y1} must be
  numeric vectors of equal length. They are interpreted as the cartesian
  coordinates of the endpoints of the line segments.

  A line segment pattern is assumed to have been observed within a specific
  region of the plane called the observation window.
  An object of class \code{"psp"} representing a point pattern
  contains information specifying the observation window.
  This window must always be specified when creating a point pattern dataset;
  there is intentionally no default action of ``guessing'' the window
  dimensions from the data points alone. 

  The argument \code{window} must be an object of class
  \code{"owin"}. It is a full description of the window geometry,
  and could have been obtained from \code{\link{owin}} or
  \code{\link{as.owin}}, or by just extracting the observation window
  of another dataset, or by manipulating such windows.
  See \code{\link{owin}} or the Examples below.

  The optional argument \code{marks} is given if the line segment pattern
  is marked, i.e. if each line segment carries additional information.
  For example, line segments which are classified into two or more different
  types, or colours, may be regarded as having a mark which identifies
  which colour they are. 

  In the current implementation, \code{marks} must be a vector, of
  the same length as \code{x0}, which is interpreted so
  that \code{marks[i]} is the mark attached to the \eqn{i}th line segment.
  If the mark is a real number then \code{marks}
  should be a numeric vector, while if the mark takes only a finite
  number of possible values (e.g. colours or types) then
  \code{marks} should be a \code{factor}.
  
  See \code{\link{psp.object}} for a description of the
  class \code{"psp"}.

  Users would normally invoke \code{psp} to create a line segment pattern,
  and the function \code{\link{as.psp}} to convert data in another
  format into a line segment pattern.
}
\seealso{
  \code{\link{psp.object}},
  \code{\link{as.psp}},
  \code{\link{owin.object}},
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  X <- psp(runif(20), runif(20), runif(20), runif(20),  window=owin())
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
