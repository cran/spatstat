\name{spokes}
\alias{spokes}
\title{Spokes pattern of dummy points}
\description{
  Generates a pattern of dummy points in a window,
  given a data point pattern.
  The dummy points lie on the radii of circles emanating from each
  data point.
}
\usage{
 spokes(x, y, nrad = 3, nper = 3, fctr = 1.5, Mdefault = 1)
}
\arguments{
  \item{x}{
    Vector of \eqn{x} coordinates of data points.
  }
  \item{y}{
    Vector of \eqn{y} coordinates of data points.
  }
  \item{nrad}{
    Number of radii emanating from each data point.
  }
  \item{nper}{
    Number of dummy points per radius.
  }
  \item{fctr}{
    Scale factor.
    Length of largest spoke radius is \code{fctr * M}
    where \code{M} is the mean nearest neighbour distance
    for the data points.
  }
  \item{Mdefault}{
    Value of \code{M} to be used if \code{x} has length 1.
  }
}
\value{
  A list with two components \code{x} and \code{y}, which are numeric
  vectors giving the coordinates of the dummy points.
}
\details{
  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}).

  Given the point pattern specified by the vectors
  \code{x} and \code{y}, the function 
  creates a pattern of \code{nrad * nper * length(x)} dummy points.

  Around each data point \code{(x[i],y[i])} there are
  \code{nrad * nper} dummy points, lying on \code{nrad} radii
  emanating from \code{(x[i],y[i])}, with \code{nper} dummy points
  equally spaced along each radius.

  The (equal) spacing of dummy points along each radius is
  controlled by the factor \code{fctr}. 
  The distance from a data point to the furthest of its associated
  dummy points is \code{fctr * M}
  where \code{M} is the mean nearest neighbour distance
  for the data points.

  If the point pattern specified by \code{x} and \code{y} contains
  only one point, the nearest neighbour distance is infinite,
  so the value \code{Mdefault} will be used in place of \code{M}.
}
\seealso{
  \code{\link{quad.object}},
  \code{\link{quadscheme}},
  \code{\link{inside.owin}},
  \code{\link{gridcentres}},
  \code{\link{stratrand}}
}
\examples{
  dat <- runifrect(10)
  \dontrun{
  plot(dat)
  }
  dum <- spokes(dat$x, dat$y)
  \dontrun{
  points(dum$x, dum$y, pch=".")
  }
  Q <- quadscheme(dat, dum)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
