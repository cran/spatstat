\name{default.dummy}
\alias{default.dummy}
\title{Generate a Default Pattern of Dummy Points}
\description{
  Generates a default pattern of dummy points
  for use in a quadrature scheme.
}
\usage{
 default.dummy(X, nx, ny, random=FALSE)
}
\arguments{
  \item{X}{
    The observed data point pattern.
    An object of class \code{"ppp"}
    or in a format recognised by \code{\link{as.ppp}()}
  }
  \item{nx, ny}{
    Integers specifying an \code{nx} by \code{ny}
    rectangular array of dummy points. 
  }
  \item{random}{
    Logical value. If \code{TRUE}, the dummy points are
    randomised.
  }
} 
\value{
  A point pattern (an object of class \code{"ppp"},
  see \code{\link{ppp.object}}) containing the dummy points.
}
\details{
  This function provides a sensible default for the dummy points
  in a quadrature scheme.

  A quadrature scheme consists of 
  the original data point pattern, an additional pattern of dummy points,
  and a vector of quadrature weights for all these points.
  See \code{\link{quad.object}} for further information about
  quadrature schemes.

  If \code{random} is false (the default),
  then the function creates dummy points
  in an \code{nx} by \code{ny} rectangular grid.
  If \code{random} is true, then the frame of the window is divided into
  an  \code{nx} by \code{ny} array of tiles, and one dummy point
  is generated at random inside each tile.
  In either case, the four corner points of the frame of the window
  are added. Then if the window is not rectangular, any dummy points
  lying outside it are deleted. 
  
  If \code{nx} is missing, a default value (depending on the
  data pattern \code{X}) is computed by \code{default.ngrid}.
  If \code{ny} is missing, it defaults to \code{nx}.
  
  Alternative functions for creating dummy patterns
  include \code{\link{corners}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}} and
  \code{\link{spokes}}.
 
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}},
  \code{\link{corners}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}},
  \code{\link{spokes}}
}
\examples{
  \dontrun{
  library(spatstat)
  }
  data(simdat)
  P <- simdat
  D <- default.dummy(P, 100)
  \dontrun{plot(D)}
  Q <- quadscheme(P, D, "grid")
  \dontrun{plot(union.quad(Q))}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
