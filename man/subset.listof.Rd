\name{subset.listof}
\alias{[.listof}
\alias{[<-.listof}
\title{Extract or Replace Subset of a List of Things}
\description{
  Extract or replace a subset of a list of things.
}
\usage{
  \method{[}{listof}(x, i)
  \method{[}{listof}(x, i) <- value
}
\arguments{
  \item{x}{
    An object of class \code{"listof"} representing a list of things
    which all belong to one class.
  }
  \item{i}{
    Subset index. Any valid subset index in the usual \R sense.
  }
  \item{value}{
    Replacement value for the subset. 
  }
}
\value{
  Another object of class \code{"listof"}.
}
\details{
  These are subset methods for the class \code{"listof"}.

  The argument \code{x} should be an object of class \code{"listof"}
  representing a list of things that all belong to one class.

  The methods extract or replace a designated
  subset of \code{x}, and return an object of class \code{"listof"}.
}
\seealso{
  \code{\link{plot.listof}},
  \code{\link{summary.listof}}
}
\examples{
   x <- list(A=runif(10), B=runif(10), C=runif(10))
   class(x) <- c("listof", class(x))
   x[1:2]
   x[1] <- list(A=rnorm(10))
 }
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{manip}
