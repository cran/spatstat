\name{pairdist.default}
\alias{pairdist.default}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a set of points
}
\usage{
  \method{pairdist}{default}(X, Y=NULL, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    Arguments specifying the coordinates of a set of points.
    Typically \code{X} and \code{Y} would be
    numeric vectors of equal length.
    Alternatively \code{Y} may be omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored.
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
  Given the coordinates of a set of points,
  this function computes the Euclidean distances between all pairs of
  points, and returns the matrix of distances.
  It is a method for the generic function \code{pairdist}.

  The arguments \code{X} and \code{Y} must determine
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
   x <- runif(100)
   y <- runif(100)
   d <- pairdist(x, y)
   d <- pairdist(cbind(x,y))
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and 
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}
\keyword{math}

