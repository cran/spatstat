\name{fv.object}
\alias{fv.object}
\title{Data Frames of Function Values}
\description{
  A class \code{"fv"} to support the convenient plotting
  of several estimates of the same function.
}
\details{
  An object of this class is a convenient way of storing and plotting
  several different estimates of the same function.

  It is a data frame with extra attributes indicating
  the recommended way of plotting the function, and other information.

  There are methods for \code{print} and \code{plot} for
  this class.
  
  Objects of class \code{"fv"} are returned by
  \code{\link{Fest}}, \code{\link{Gest}},\code{\link{Jest}},
  and \code{\link{Kest}}.
}
\seealso{
  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Jest}},
  \code{\link{Kest}}
}
\examples{
    require(spatstat)
    data(cells)
    K <- Kest(cells)

    class(K)

    K  # prints a sensible summary

    plot(K)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
 
