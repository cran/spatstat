\name{rsyst}
\alias{rsyst}
\title{Systematic random point pattern}
\description{
  Generates a ``systematic random'' pattern of points in a window,
  consisting of a grid of equally-spaced points with a random common
  displacement.
}
\usage{
 rsyst(win=square(1), nx, ny=nx, dx=NULL, dy=NULL)
}
\arguments{
  \item{win}{
    A window. 
    An object of class \code{\link{owin}},
    or data in any format acceptable to \code{\link{as.owin}()}.
  }
  \item{dx}{Spacing of grid points in \eqn{x} direction.
    Incompatible with \code{nx}.
  }
  \item{dy}{Spacing of grid points in \eqn{y} direction.
    Incompatible with \code{ny}.
  }
  \item{nx}{Number of columns of grid points in the window.
    Incompatible with \code{dx}.
  }
  \item{ny}{Number of rows of grid points in the window.
    Incompatible with \code{dy}.
  }
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\details{
  This function generates a \dQuote{systematic random} pattern of points
  in the window \code{win}. The pattern consists of a rectangular grid
  of points with a random common displacement.

  The grid spacing is determined by the distances \code{dx, dy}
  if they are present. If they are absent, then the grid spacing is
  determined so that there will be \code{nx} columns and \code{ny} rows
  of grid points in the bounding rectangle of \code{win}.

  The grid is then given a random displacement (the common displacement
  of the grid points is a uniformly distributed random vector in the
  tile of dimensions \code{dx, dy}).

  Some of the resulting grid points may lie outside the window \code{win}:
  if they do, they are deleted.
  The result is a point pattern inside the window \code{win}.

  This function is useful in creating dummy points for quadrature
  schemes (see \code{\link{quadscheme}}) as well as in simulating
  random point patterns.
}
\seealso{
  \code{\link{rstrat}},
  \code{\link{runifpoint}},
  \code{\link{quadscheme}}
}
\examples{
  X <- rsyst(nx=10)
  plot(X)

  # polygonal boundary
  data(letterR)
  X <- rsyst(letterR, 5, 10)
  plot(X)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{datagen}

