\name{pairdist}
\alias{pairdist}
\alias{pairdist.ppp}
\alias{pairdist.default}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a point pattern.
}
\usage{
  pairdist(X, \dots, method="C")
  pairdist.ppp(X, \dots, method="C")
  pairdist.default(X, Y=NULL, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    Arguments specifying the locations of a set of points.
    For \code{pairdist.ppp}, the argument \code{X} should be a point
    pattern (object of class \code{"ppp"}).
    For \code{pairdist.default}, typically \code{X} and \code{Y} would be
    numeric vectors of equal length. Alternatively \code{Y} may be
    omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored by \code{pairdist.ppp}
    and \code{pairdist.default}.
  }
  \item{method}{
    String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
    Usually not specified.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the points numbered \code{i} and \code{j}.
}
\details{
  This function computes the Euclidean distances between all pairs of
  points, and returns the matrix of distances.

  The function \code{pairdist} is generic, with
  a method for point patterns (objects of class \code{"ppp"})
  and a default method.

  The method for point patterns expects a single
  point pattern argument \code{X} and returns the matrix of its
  pairwise distances.

  The default method expects that \code{X} and \code{Y} will determine
  the coordinates of a set of points. Typically \code{X} and
  \code{Y} would be numeric vectors of equal length. Alternatively
  \code{Y} may be omitted and \code{X} may be a list with two components
  named \code{x} and \code{y}, or a matrix or data frame with two columns.
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. The C code is somewhat faster.
}
\seealso{
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\examples{
   data(cells)
   d <- pairdist(cells)

   x <- runif(100)
   y <- runif(100)
   d <- pairdist(x, y)
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and 
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

