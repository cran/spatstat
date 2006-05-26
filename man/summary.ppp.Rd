\name{summary.ppp}
\alias{summary.ppp}
\title{Summary of a Point Pattern Dataset}
\description{
  Prints a useful summary of a point pattern dataset.
}
\usage{
 \method{summary}{ppp}(object, \dots, checkdup=TRUE)
}
\arguments{
  \item{object}{
    Point pattern (object of class \code{"ppp"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{checkdup}{
    Logical value indicating whether to check
    for the presence of duplicate points.
  }
}
\details{
  A useful summary of the point pattern \code{object} is printed.

  This is a method for the generic function \code{\link{summary}}.

  If \code{checkdup=TRUE}, the pattern will be checked for the
  presence of dublicate points, using \code{\link{duplicated.ppp}}.
  This can be time-consuming if the pattern contains many points,
  so the checking can be disabled by setting \code{checkdup=FALSE}.
}
\seealso{
  \code{\link{summary}},
  \code{\link{summary.owin}},
  \code{\link{print.ppp}}
}
\examples{
  data(cells)      # plain vanilla point pattern
  summary(cells)  

  data(lansing)    # multitype point pattern
  \testonly{lansing <- lansing[seq(1, lansing$n, length=40)]}
  summary(lansing) # tabulates frequencies of each mark
  
  data(longleaf)    # numeric marks
  \testonly{longleaf <- longleaf[seq(1, longleaf$n, length=40)]}
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
\keyword{methods}

