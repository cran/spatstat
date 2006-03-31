\name{nndist}
\alias{nndist}
\alias{nndist.ppp}
\alias{nndist.default}
\title{Nearest neighbour distances}
\description{
  Computes the distance from each point to its nearest neighbour
  in a point pattern.
}
\usage{
  nndist(X, \dots, method="C")
  nndist.ppp(X, \dots, method="C")
  nndist.default(X, Y=NULL, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    Arguments specifying the locations of
    a set of points.
    For \code{nndist.ppp}, the argument \code{X} should be a point
    pattern (object of class \code{"ppp"}).
    For \code{nndist.default}, typically \code{X} and \code{Y} would be
    numeric vectors of equal length. Alternatively \code{Y} may be
    omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored by \code{nndist.ppp}
    and \code{nndist.default}.
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

  The function \code{nndist} is generic, with
  a method for point patterns (objects of class \code{"ppp"})
  and a default method.

  The method for point patterns expects a single
  point pattern argument \code{X} and returns the vector of its
  nearest neighbour distances.

  The default method expects that \code{X} and \code{Y} will determine
  the coordinates of a set of points. Typically \code{X} and
  \code{Y} would be numeric vectors of equal length. Alternatively
  \code{Y} may be omitted and \code{X} may be a list with two components
  named \code{x} and \code{y}, or a matrix or data frame with two columns.
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="C"}
  (the default) then C code is used. 
  The C code is faster by two to three orders of magnitude
  and uses much less memory.
  
  If there is only one point (if \code{x} has length 1),
  then a nearest neighbour distance of \code{Inf} is returned.
  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.
}
\section{Warnings}{
  An infinite value is returned if there is only one point
  in the point pattern.
}
\seealso{
  \code{\link{pairdist}},
  \code{\link{Gest}}
}
\examples{
   data(cells)
   d <- nndist(cells)

   x <- runif(100)
   y <- runif(100)
   d <- nndist(x, y)

   # Stienen diagram
   plot(cells \%mark\% (nndist(cells)/2), markscale=1)
}
\author{Pavel Grabarnik
  \email{pavel.grabar@issp.serpukhov.su}
  and
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
}
\keyword{spatial}

