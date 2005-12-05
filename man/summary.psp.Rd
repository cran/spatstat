\name{summary.psp}
\alias{summary.psp}
\title{Summary of a Line Segment Pattern Dataset}
\description{
  Prints a useful summary of a line segment pattern dataset.
}
\synopsis{
  summary.psp(object, \dots)
}
\usage{
 summary.psp(object)
}
\arguments{
  \item{object}{Line segment pattern (object of class \code{"psp"}).}
}
\details{
  A useful summary of the line segment pattern \code{object} is printed.

  This is a method for the generic function \code{\link{summary}}.
}
\seealso{
  \code{\link{summary}},
  \code{\link{summary.owin}},
  \code{\link{print.psp}}
}
\examples{
  a <- psp(runif(10), runif(10), runif(10), runif(10), window=owin())
  summary(a)  # describes it
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
