\name{nndist}
\alias{nndist}
\title{Nearest neighbour distances}
\description{
  Computes the distance from each point to its nearest neighbour
  in a point pattern.
}
\usage{
  nndist(x, y=NULL, method="C")
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
  Numeric vector of the nearest neighbour distances for each point.
}
\details{
  This function computes the Euclidean distance from each point
  in a point pattern to its nearest neighbour (the nearest other
  point of the pattern).

  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  then C code is used. The C code is faster by two to three orders of magnitude
  and uses much less memory.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{Gest}}
}
\examples{
   x <- runif(100)
   y <- runif(100)
   d <- nndist(x, y)
   data(cells)
   d <- nndist(cells)
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

