\name{nndist.psp}
\alias{nndist.psp}
\title{Nearest neighbour distances between line segments}
\description{
  Computes the distance from each line segment 
  to its nearest neighbour in a line segment pattern.
  Alternatively finds the distance to the second nearest,
  third nearest etc.
}
\usage{
  \method{nndist}{psp}(X, \dots, k=1, method="Fortran")
}
\arguments{
  \item{X}{
    A line segment pattern (object of class \code{"psp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{k}{
    Integer. The algorithm will find
    the distance to the \code{k}th nearest neighbour.
  }
  \item{method}{
    String specifying which method of calculation to use.
    Values are \code{"Fortran"} and \code{"interpreted"}.
    Usually not specified.
  }
}
\value{
  A numeric vector whose \code{[i]} entry is the distance
  from the line segment numbered \code{i} to the nearest other line
  segment (or in general the distance to the \code{k}th nearest line segment).
}
\details{
  If \code{k=1}, this function computes the distance from each line segment 
  to the nearest other line segment in \code{X}.
  In general it computes the distance from each line segment to the
  \code{k}th nearest other line segment.

  This is a method for the generic function \code{\link{nndist}}
  for the class \code{"psp"}.

  Distances are calculated using the Hausdorff metric. The Hausdorff
  distance between two line segments is the maximum distance
  from any point on one of the segments to the nearest point on
  the other segment.   

  If there are fewer than \code{k} line segments in the pattern,
  the \code{k}th nearest neighbour distances will be infinite (\code{Inf}).
  
  The argument \code{method} is not normally used. It is
  retained only for checking the validity of the software.
  If \code{method = "interpreted"} then the distances are
  computed using interpreted R code only. If \code{method="Fortran"}
  (the default) then Fortran code is used. The Fortran code is somewhat faster.
}
\seealso{
  \code{\link{nndist}},
  \code{\link{nndist.ppp}}
}
\examples{
   L <- psp(runif(10), runif(10), runif(10), runif(10), owin())
   D <- nndist(L)
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}

