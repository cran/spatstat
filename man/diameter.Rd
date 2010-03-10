\name{diameter}
\alias{diameter}
\title{Diameter of a Window}
\description{
  Computes the diameter of a window 
}
\usage{
 diameter(w)
}
\arguments{
  \item{w}{A window, whose diameter will be computed.}
}
\value{
  The numerical value of the diameter of the window. 
}
\details{
  This function computes the 
  diameter of a window of arbitrary shape,
  i.e. the maximum distance between any two points
  in the window. 

  The argument \code{w} should be a window (an object of class
  \code{"owin"}, see \code{\link{owin.object}} for details)
  or can be given in any format acceptable to \code{\link{as.owin}()}.
}
\seealso{
  \code{\link{area.owin}},
  \code{\link{perimeter}},
  \code{\link{owin}},
  \code{\link{as.owin}}
}
\examples{
  w <- owin(c(0,1),c(0,1))
  diameter(w) 
  # returns sqrt(2)
  data(letterR)
  diameter(letterR)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{math}
