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
    an object of class \code{"spp"} as defined in the 
    \code{spatial} library
    \item
    an object of class \code{"quad"} representing a quadrature scheme
    (see \code{\link{quad.object}})
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

  An \code{spp} object is the representation of a point pattern
  in the \code{spatial} library. Our implementation recognises
  the following formats:
  \itemize{
    \item
    a structure with entries \code{x}, \code{y}
    \code{xl}, \code{xu}, \code{yl}, \code{yu}
    \item
    a structure with entries
    \code{x}, \code{y} and \code{area},
    where \code{area} is a structure with entries
    \code{xl}, \code{xu}, \code{yl}, \code{yu}
  }
  (used in \code{spatial} versions 1 to 6 and version 7.1 respectively)
  where \code{x} and \code{y} are vectors of equal length
  giving the point coordinates, and \code{xl}, \code{xu}, \code{yl},
  \code{yu} are numbers giving the dimensions of a rectangular window.

  Point pattern datasets can also be created by the function
  \code{\link{ppp}}.
}
\seealso{
  \code{\link{ppp}},
  \code{\link{ppp.object}},
  \code{\link{as.owin}},
  \code{\link{owin.object}}
}
\section{Warnings}{
  If the format of \code{spp} objects
  is changed in future versions of the \code{spatial} library,
  then \code{as.ppp} may not be able to interpret them.
  It currently handles all versions of \code{spatial} up to 7.1-4.
}
\examples{
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
 require(spatial)
 towns <- ppinit("towns.dat")
 pp <- as.ppp(towns)   # converted to our format
 detach(package:spatial)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
