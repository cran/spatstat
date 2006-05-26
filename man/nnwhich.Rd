\name{nnwhich}
\alias{nnwhich}
\alias{nnwhich.ppp}
\alias{nnwhich.default}
\title{Nearest neighbour}
\description{
  Finds the nearest neighbour of each point in a point pattern.
}
\usage{
  nnwhich(X, \dots, method="C")
  \method{nnwhich}{ppp}(X, \dots, method="C")
  \method{nnwhich}{default}(X, Y=NULL, \dots, method="C")
}
\arguments{
  \item{X,Y}{
    Arguments specifying the locations of
    a set of points.
    For \code{nnwhich.ppp}, the argument \code{X} should be a point
    pattern (object of class \code{"ppp"}).
    For \code{nnwhich.default}, typically \code{X} and \code{Y} would be
    numeric vectors of equal length. Alternatively \code{Y} may be
    omitted and \code{X} may be
    a list with two components \code{x} and \code{y},
    or a matrix with two columns.
  }
  \item{\dots}{
    Ignored by \code{nnwhich.ppp}
    and \code{nnwhich.default}.
  }
  \item{method}{String specifying which method of calculation to use.
    Values are \code{"C"} and \code{"interpreted"}.
  }
}
\value{
  Integer vector giving, for each point, the index of its nearest neighour.
}
\details{
  For each point in the given point pattern, this function finds
  its nearest neighbour (the nearest other point of the pattern).
  It returns a vector giving, for each point, the index of the point's
  nearest neghbour.

  The function \code{nnwhich} is generic, with
  a method for point patterns (objects of class \code{"ppp"})
  and a default method.

  The method for point patterns expects a single
  point pattern argument \code{X}.
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
  then the nearest neighbour is undefined, and a value of \code{NA}
  is returned.
  If there are no points (if \code{x} has length zero)
  a numeric vector of length zero is returned.

  To evaluate the \emph{distance} between a point and its nearest
  neighbour, use \code{\link{nndist}}.
}
\section{Warnings}{
  A value of \code{NA} is returned if there is only one point
  in the point pattern.
}
\seealso{
  \code{\link{nndist}}
}
\examples{
   oldpar <- par(mfrow=c(2,1))
   data(cells)
   plot(cells)
   m <- nnwhich(cells)

   # plot nearest neighbour links
   b <- cells[m]
   arrows(cells$x, cells$y, b$x, b$y, angle=15, length=0.15, col="red")

   # find points which are the neighbour of their neighbour
   self <- (m[m] == seq(m))
   # plot them
   A <- cells[self]
   B <- cells[m[self]]
   plot(cells)
   segments(A$x, A$y, B$x, B$y)
   par(oldpar)
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

