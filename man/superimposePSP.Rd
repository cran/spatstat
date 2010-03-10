\name{superimposePSP}
\alias{superimposePSP}
\title{Superimpose Several Line Segment Patterns}
\description{
  Superimpose any number of line segment patterns.
}
\usage{
  superimposePSP(\dots, W=NULL, check=TRUE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which is a line segment pattern
    (object of class \code{"psp"}).
  }
  \item{W}{
    Optional. Window for the resulting line segment pattern.
    An object of class \code{"owin"}, or something
    acceptable to \code{\link{as.owin}}.
  }
  \item{check}{
    Logical value (passed to \code{\link{psp}})
    determining whether to check the geometrical validity of the
    resulting line segment pattern.
  }
}
\value{
  A line segment pattern (object of class \code{"psp"}) representing the
  superposition (union) of all the line segment patterns.
}
\details{
  This function is used to superimpose two or more line segment patterns.
  
  Each of the arguments in
  \code{\dots} is a line segment pattern object (of class \code{"psp"}).
  The line segment patterns are not required to have
  the same window of observation.

  The window for the superimposed line segment pattern is
  specified by \code{W}. Its default value is
  the union of the windows of all the line segment patterns.

  If any of the arguments is a \emph{marked} line segment pattern,
  then all of them must be.
  In that case, the result is also a marked line segment pattern.
}
\seealso{
  \code{\link{superimpose}} for superimposing point patterns.
}
\examples{
   X <- rpoisline(10)
   Y <- as.psp(matrix(runif(40), 10, 4), window=owin())
   Z <- superimposePSP(X, Y)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
