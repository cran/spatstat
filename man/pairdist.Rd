\name{pairdist}
\alias{pairdist}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a point pattern.
}
\usage{
  pairdist(x, y=NULL, method="C")
}
\arguments{
  \item{x,y}{
    Arguments specifying the \eqn{x} and \eqn{y} coordinates of
    a set of points. Typically \code{x} and \code{y} would be
    numeric vectors of equal length. Alternatively \code{y} may be
    omitted and \code{x} may be
    a point pattern (object of class \code{"ppp"})
    or a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  A square matrix whose \code{[i,j]} entry is the distance
  between the points numbered \code{i} and \code{j}.
}
\details{
  This function computes the Euclidean distances between all pairs of
  points, and returns the matrix of distances.

  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  then C code is used. The C code is marginally faster.
}
\seealso{
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\examples{
   x <- runif(100)
   y <- runif(100)
   d <- pairdist(x, y)
   data(cells)
   d <- pairdist(cells)
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and 
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

