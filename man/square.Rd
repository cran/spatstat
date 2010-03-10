\name{square}
\alias{square}
\alias{unit.square}
\title{Square Window}
\description{
Creates a square window
}
\usage{
 square(r=1)
 unit.square()
}
\arguments{
  \item{r}{The side length of the square.}
}
\value{
  An object of class \code{"owin"} (see \code{\link{owin.object}})
  specifying a window. 
}
\details{
  \code{square}
  is a shortcut for creating a window object
  representing the square 
  \eqn{[0,r] \times [0,r]}{[0,r] * [0,r]}.
  It is equivalent to the command
  \code{owin(c(0,r),c(0,r))}.

  \code{unit.square} creates the unit square
  \eqn{[0,1] \times [0,1]}{[0,1] * [0,1]}.
  It is equivalent to
  \code{square(1)} or \code{square()} or \code{owin(c(0,1),c(0,1))}.

  These commands are 
  included mainly to improve the readability of some code.
}
\seealso{
  \code{\link{owin.object}},
  \code{\link{owin}}
}
\examples{
 W <- square(10)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
 
 
