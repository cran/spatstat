\name{unmark}
\alias{unmark}
\title{Remove Marks from a Marked Point Pattern}
\description{
  Removes the mark information from a marked point pattern
}
\usage{
 unmark(X)
}
\arguments{
  \item{X}{A point pattern. Either an object of class \code{"ppp"}
    or anything acceptable to \code{\link{as.ppp}}.
  }
}
\value{
  A point pattern (of class \code{"ppp"})
  with any mark information deleted.
}
\details{
  This is simply a convenience for those occasions in which we want to
  remove the marks from a marked point pattern.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{as.ppp}}
}
\examples{
  \dontrun{library(spatstat)}
  data(lansing)
  hicks <- lansing[lansing$marks == "hickory", ]
  \dontrun{
  plot(hicks)  # still a marked point pattern, but only 1 value of marks
  plot(unmark(hicks)) # unmarked
  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
