\name{summary.quad}
\alias{summary.quad}
\alias{print.summary.quad}
\title{Summarizing a Quadrature Scheme}
\description{
  \code{summary} method for class \code{"quad"}.
}
\usage{
  summary.quad(object, \dots)
  print.summary.quad(x, \dots, dp=3)
}
\arguments{
  \item{object}{A quadrature scheme.}
  \item{\dots}{Ignored.}
  \item{dp}{Number of significant digits to print.} 
  \item{x}{Object of class \code{"summary.quad"} returned by
    \code{summary.quad}.}
}
\details{
  This is a method for the generic \code{\link{summary}}
  for the class \code{"quad"}. An object of class \code{"quad"}
  describes a quadrature scheme, used to fit a point process model.
  See \code{\link{quad.object}}) for details of this class. 
  
  \code{summary.quad} extracts information about the
  quadrature scheme,
  and \code{print.summary.quad} prints this information in a
  comprehensible format.

  In normal usage, \code{print.summary.quad} is invoked implicitly
  when the user calls \code{summary.quad} without assigning its value
  to anything. See the examples.
}
\value{
  \code{summary.quad} returns an object of class \code{"summary.quad"},
  while \code{print.summary.quad} returns \code{NULL}.
}
\examples{
  require(spatstat)
  # make a quadrature scheme
  Q <- quadscheme(rpoispp(42))
  # summarize it
  summary(Q)
  # save the summary
  s <- summary(Q)
  # print it
  print(s)
  s
  # extract total quadrature weight
  s$w$all$sum
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
