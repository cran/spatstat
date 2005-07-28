\name{compatible.fv}
\alias{compatible.fv}
\title{Test Whether Two Function Objects Are Compatible}
\description{
  Tests whether two function objects (class \code{"fv"})
  are compatible.
}
\usage{
  compatible.fv(A, B)
}
\arguments{
  \item{A,B}{Function value objects (class \code{"fv"}).}
}
\details{
  An object of class \code{"fv"} is essentially a data frame
  containing several different statistical estimates of the same
  function. Such objects are returned by \code{\link{Kest}} and its
  relatives.

  This function tests whether two such objects are compatible
  (so that, for example, they could be added or subtracted).
  
  The functions are compatible if they have been evaluated at the
  same sequence of values of the argument \code{r}, and if the
  statistical estimates have the same names.
}
\value{
  Logical value: \code{TRUE} if the objects are compatible, and \code{FALSE}
  if they are not.
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
