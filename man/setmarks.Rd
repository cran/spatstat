\name{setmarks}
\alias{setmarks}
\alias{"%mark%"}
\title{Set or Reset the Marks in a Point Pattern}
\description{
  Create a marked point pattern from given point locations
  and given marks.
}
\synopsis{
  setmarks(X, m)
}
\usage{
  setmarks(X, m)
  X %mark% m
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"})
  }
  \item{m}{
    Vector of mark values (of any atomic mode)
  }
}
\value{
  Marked point pattern (object of class \code{"ppp"})
  identical to \code{X} except that it has marks equal to \code{m}.
}
\details{
  This function is a shortcut which creates a marked point pattern
  from a given point pattern \code{X} and a given vector \code{m} 
  by attaching the marks in \code{m} to the locations of the points
  in \code{X}.

  If \code{X} is already a marked point pattern, then its marks
  are ignored and replaced by the values in \code{m}.
  
  The length of \code{m} must agree with the number of points in
  \code{X}. 

  Use \code{\link{unmark}} to remove marks.
  Use \code{\link{ppp}} to create point patterns in more general
  situations.
}
\seealso{
  \code{\link{unmark}},
  \code{\link{ppp.object}},
  \code{\link{is.marked}},
  \code{\link{is.marked.ppp}}
}
\examples{
   require(spatstat)
   data(cells)

   m <- runif(cells$n)

   Y <- setmarks(cells, m)
   Y <- cells %mark% m
   # equivalent

   \dontrun{plot(Y)}
   is.marked(Y)  #TRUE
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
