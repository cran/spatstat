\name{summary.ppp}
\alias{summary.ppp}
\title{Summary of a Point Pattern Dataset}
\description{
  Prints a useful summary of a point pattern dataset.
}
\synopsis{
  summary.ppp(object, \dots)
}
\usage{
 summary.ppp(object)
 summary(object)
}
\arguments{
  \item{object}{Point pattern (object of class \code{"ppp"}).}
}
\details{
  A useful summary of the point pattern \code{object} is printed.

  This is a method for the generic function \code{\link{summary}}.
}
\seealso{
  \code{\link{summary}},
  \code{\link{summary.owin}},
  \code{\link{print.ppp}}
}
\examples{
  \dontrun{library(spatstat)}

  data(cells)      # plain vanilla point pattern
  summary(cells)  

  data(lansing)    # multitype point pattern
  summary(lansing) # tabulates frequences of each mark
  
  data(longleaf)    # numeric marks
  summary(longleaf) # prints summary.default(x$marks)

  data(demopat)     # weird polygonal window
  summary(demopat)  # describes it
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
