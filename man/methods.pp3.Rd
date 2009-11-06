\name{methods.pp3}
\Rdversion{1.1}
\alias{print.pp3}
\alias{unitname.pp3}
\title{
  Methods for three-dimensional point patterns
}
\description{
  Methods for class \code{"pp3"}.
}
\usage{
  \method{print}{pp3}(x, ...)
  \method{unitname}{pp3}(x)
}
\arguments{
  \item{x}{
    Object of class \code{"pp3"}.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  These are methods for the generic functions
  \code{\link{print}} and \code{\link{unitname}}
  for the class \code{"pp3"} of three-dimensional point patterns.

  The \code{print} method prints a description of the point pattern,
  while the \code{unitname} method extracts the name of the unit of
  length in which the point coordinates are expressed.
}
\value{
  For \code{print.pp3} the value is \code{NULL}.
  For \code{unitname.pp3} an object of class \code{"units"}.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{pp3}},
  \code{\link{print}},
  \code{\link{unitname}}
}
\examples{
   X <- pp3(runif(42),runif(42),runif(42), box3(c(0,1), unitname="mm"))
   X
   unitname(X) 
}
\keyword{spatial}
\keyword{methods}

