\name{as.ppp}
\alias{as.ppp}
\title{Convert Data To Class ppp}
\description{
Tries to coerce any reasonable kind of data to a point pattern
(an object of class \code{"ppp"})
for use by the \code{spatstat} package).
}
\synopsis{
  as.ppp(X, W=NULL, fatal=TRUE)
}
\usage{
  as.ppp(X)
  as.ppp(X, W)
  as.ppp(X, fatal=TRUE)
}
\arguments{
  \item{X}{Data which will be converted into a point pattern}
  \item{W}{Data which define a default window for the pattern}
  \item{fatal}{Logical value. See Details.}
}
\value{
  An object of class \code{"ppp"} (see \code{\link{ppp.object}})
  describing the point pattern and its window of observation.
  The value \code{NULL} may also be returned; see Details.
}
\details{
  This function tries to coerce the data \code{X} to a point pattern
  (an object of class \code{"ppp"}; see \code{\link{ppp.object}} for
  an overview).  The optional argument \code{W} is a default window for the
  pattern; it is used only if a window is not specified in \code{X}.

  The data \code{X} may be:
  \itemize{
    \item
    an object of class \code{"ppp"} 
    \item
    an object of class \code{"spp"} as used in the 
    Venables and Ripley package \code{spatial}
    \item
    an object of class \code{"quad"} representing a quadrature scheme
    (see \code{\link{quad.object}})
    \item
    a structure with entries \code{x}, \code{y} (vectors of equal
    length) and \code{xl}, \code{xu}, \code{yl}, \code{yu} (numbers
    giving the dimensions of a rectangular window)
    \item
    a two-column matrix
    \item
    a structure with entries \code{x}, \code{y}
  }
  In the last two cases, we need the second argument \code{W}
  which is converted to a window object
  by the function \code{\link{as.owin}}.
  In the first three cases, \code{W} will be ignored.

  The argument \code{fatal} indicates what to do when 
  \code{W} is missing and \code{X} contains no
  information about the window. If \code{fatal=TRUE}, a fatal error
  will be generated; if \code{fatal=FALSE}, the
  value \code{NULL} is returned.

  Point pattern datasets can also be created by the function
  \code{\link{ppp}}.
}
\seealso{
  \code{\link{ppp}},
  \code{\link{ppp.object}},
  \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\examples{
 library(spatstat)

 \dontrun{
 plot(c(0,1),c(0,1),type="n")
 xy <- locator(20) # click on 20 points in plot window
 pp <- as.ppp(xy, c(0,1,0,1))
 
 w <- owin(c(0,1),c(0,1))
 plot(w)                      # neater
 xy <- locator(20)            # click on 20 points in plot window
 pp <- as.ppp(xy, w)
 } 
 
 xy <- matrix(runif(40), ncol=2)
 pp <- as.ppp(xy, c(0,1,0,1))
 
 # Venables-Ripley format
 sp <- list(x=runif(10), y=runif(10), xl=0, xu=1, yl=0, yu=1)
 pp <- as.ppp(sp)   # converted to our format
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
