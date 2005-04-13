\name{crossdist}
\alias{crossdist}
\alias{crossdist.ppp}
\alias{crossdist.default}
\title{Pairwise distances between two different point patterns}
\description{
  Computes the distances between pairs of points
  taken from two different point patterns.
}
\usage{
  crossdist(X, Y, \dots, method="C")
  crossdist.ppp(X, Y, \dots, method="C")
  crossdist.default(X, Y, x2, y2, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    For \code{crossdist.ppp}, the arguments \code{X} and \code{Y}
    should be point patterns (objects of class \code{"ppp"}).
    For \code{crossdist.default}, the arguments \code{X} and \code{Y}
    should be numeric vectors of equal length specifying the location of
    the first point pattern.
  }
  \item{x2,y2}{
    For \code{crossdist.default}, these
    should be numeric vectors of equal length specifying the location of
    the second point pattern.
  }
  \item{\dots}{
    Ignored by \code{crossdist.ppp}
    and \code{crossdist.default}.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  A matrix whose \code{[i,j]} entry is the distance
  from the \code{i}-th point in the first pattern
  to the \code{j}-th point in the second pattern.
}
\details{
  Given two point patterns, 
  this function computes the Euclidean distance from each point
  in the first pattern to each point in the second pattern,
  and returns a matrix containing these distances.

  The function \code{crossdist} is generic, with
  a method for point patterns (objects of class \code{"ppp"})
  and a default method.

  The method for point patterns expects two
  point patterns \code{X} and \code{Y}, and returns the matrix
  whose \code{[i,j]} entry is the distance from \code{X[i]} to
  \code{Y[j]}.

  The default method expects \code{X} and \code{Y} to be numeric vectors
  of equal length specifying the coordinates of the first point pattern.
  The arguments \code{x2},\code{y2} specify the coordinates of the
  second point pattern.
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is faster by a factor of 4.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{nndist}},
  \code{\link{Gest}}
}
\examples{
   data(cells)
   d <- crossdist(cells, runifpoint(6))
   
   d <- crossdist(runif(7), runif(7), runif(12), runif(12))
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

