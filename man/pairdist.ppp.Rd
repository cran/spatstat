\name{pairdist.ppp}
\alias{pairdist.ppp}
\title{Pairwise distances}
\description{
  Computes the matrix of distances between all pairs of points
  in a point pattern.
}
\usage{
  \method{pairdist}{ppp}(X, \dots, method="C")
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
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
  This is a method for the generic function \code{pairdist}.

  Given a point pattern \code{X} (an object of class \code{"ppp"}),
  this function computes the Euclidean distances between all pairs of
  points in \code{X}, and returns the matrix of distances.

  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. The C code is somewhat faster.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{pairdist.default}},
  \code{\link{pairdist.psp}},
  \code{\link{crossdist}},
  \code{\link{nndist}},
  \code{\link{Kest}}
}
\examples{
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

