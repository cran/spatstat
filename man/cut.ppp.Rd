\name{cut.ppp}
\alias{cut.ppp}
\title{Convert Point Pattern Marks from Numeric to Factor}
\description{
  Transform the marks of a marked point pattern
  from numeric values into a factor.
}
\synopsis{
  cut.ppp(x, \dots)
}
\usage{
  cut(x, \dots)
  cut.ppp(x, \dots)
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
    This should be a marked point pattern, and the marks would normally
    be numerical values.
  }
  \item{\dots}{
    Arguments passed to the generic \code{\link{cut}} function.
    They determine the breakpoints for the mapping from numerical values to
    factor values. See \code{\link{cut}}.
  }
} 
\value{
  A multitype point pattern, that is, a point pattern object
  (of class \code{"ppp"}) with a \code{marks} vector that is a factor.
}
\details{
  This simple function applies the generic \code{\link{cut}} operation
  to the marks of the point pattern \code{x}. The range of values
  of the numerical marks is divided into several intervals, and each
  interval is associated with a level of a factor. 
  The result is a
  marked point pattern, with the same window and point locations as
  \code{x}, but with the numeric mark of each point discretised
  by replacing it by the factor level. 

  This function is a convenient
  way to transform a marked point pattern which has numeric marks
  into a multitype point pattern,
  for example to plot it or analyse it. See the examples.

  To select some points from a point pattern, use the subset operator
  \code{\link{"[.ppp"}} instead.
}
\seealso{
  \code{\link{cut}},
  \code{\link{ppp.object}}
}
\examples{
 \dontrun{library(spatstat)}

 data(longleaf)
 # Longleaf Pines data
 # the marks are positive real numbers indicating tree diameters.

 \dontrun{
 plot(longleaf)
 }

 # cut the range of tree diameters into three intervals
 long3 <- cut(longleaf, 3)
 \dontrun{
 plot(long3)
 }

 # adult trees defined to have diameter at least 30 cm
 long2 <- cut(longleaf, breaks=c(0,30,100), labels=c("Sapling", "Adult"))
 plot(long2)
 plot(long2, cols=c("green","blue"))
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
