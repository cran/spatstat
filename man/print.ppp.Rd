\name{print.ppp}
\alias{print.ppp}
\title{Print Brief Details of a Point Pattern Dataset}
\description{
  Prints a very brief description of a point pattern dataset.
}
\synopsis{
  \method{print}{ppp}(x, \dots)
}
\usage{
 \method{print}{ppp}(x)
 x
}
\arguments{
  \item{x}{Point pattern (object of class \code{"ppp"}).}
}
\details{
  A very brief description of the point pattern \code{x} is printed.

  This is a method for the generic function \code{\link{print}}.
}
\seealso{
  \code{\link{print}},
  \code{\link{print.owin}},
  \code{\link{summary.ppp}}
}
\examples{
  data(cells)      # plain vanilla point pattern
  cells

  data(lansing)    # multitype point pattern
  lansing          
  
  data(longleaf)    # numeric marks
  longleaf          

  data(demopat)     # weird polygonal window
  demopat
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{print}

