\name{as.owin}
\alias{as.owin}
\title{Convert Data To Class owin}
\description{
Converts data specifying an observation window
in any of several formats, into an object of class \code{"owin"}.
}
\usage{
 as.owin(W)
}
\arguments{
  \item{W}{Data specifying an observation window, in any of several formats
    described under \emph{Details} below.}
}
\value{
  An object of class \code{"owin"} (see \code{\link{owin.object}})
  specifying an observation window.
}
\details{
  The class \code{"owin"} is a way of specifying the observation window
  for a point pattern. See \code{\link{owin.object}} for an overview.
 
  This function converts data in any of several formats 
  into an object of class \code{"owin"} for use by the \code{spatstat}
  package.  The argument \code{W} may be
  \itemize{
    \item
    an object of class \code{"owin"}
    \item
    a structure with entries \code{xrange}, \code{yrange} specifying the 
    \eqn{x} and \eqn{y} dimensions of a rectangle
    \item
    a four-element vector
    (interpreted as \code{(xmin, xmax, ymin, ymax)})
    specifying the \eqn{x} and \eqn{y} dimensions of a rectangle
    \item
    a structure with entries \code{xl}, \code{xu}, \code{yl}, \code{yu}
    specifying the \eqn{x} and \eqn{y} dimensions of a rectangle
    as \code{(xmin, xmax) = (xl, xu)} and 
    \code{(ymin, ymax) = (yl, yu)}. This will accept objects of
    class \code{spp} used in the Venables and Ripley \code{spatial}
    library.
    \item
    an object of class \code{ppp} 
    (in which case the object's \code{window} structure will be
    extracted).
  }
}
\seealso{

  \code{\link{owin.object}},
  \code{\link{owin}}
}
\examples{
 \dontrun{library(spatstat)}
 w <- as.owin(c(0,1,0,1))
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
