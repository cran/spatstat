\name{print.owin}
\alias{print.owin}
\title{Print Brief Details of a Spatial Window}
\description{
  Prints a very brief description of a window object.
}
\synopsis{
  print.owin(x, \dots)
}
\usage{
 print.owin(x)
 print(x)
 x
}
\arguments{
  \item{x}{Window (object of class \code{"owin"}).}
}
\details{
  A very brief description of the window \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{print.ppp}},
  \code{\link{summary.owin}}
}
\examples{
  \dontrun{library(spatstat)}

  owin()  # the unit square

  data(demopat)
  W <- demopat$window  
  W                    # just says it is polygonal
  as.mask(W)           # just says it is a binary image 

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
