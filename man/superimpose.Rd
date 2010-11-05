\name{superimpose}
\alias{superimpose}
\title{Superimpose Several Point Patterns}
\description{
  Superimpose any number of point patterns.
}
\usage{
  superimpose(\dots, W=NULL, check=TRUE)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which represents a point pattern.
    Each argument must be either a point pattern
    object (of class \code{"ppp"}) or a structure
    containing elements \code{x} and \code{y}.
  }
  \item{W}{
    Optional. Window for the resulting point pattern.
    An object of class \code{"owin"}, or something
    acceptable to \code{\link{as.owin}}.
  }
  \item{check}{
    Logical value (passed to \code{\link{ppp}})
    determining whether to check the geometrical validity of the
    resulting point pattern.
  }
}
\value{
  A point pattern (object of class \code{"ppp"}) representing the
  superposition (union) of all the point patterns.
}
\details{
  This function is used to superimpose two or more point patterns
  (see also \code{\link{concatxy}}). 
  
  Each of the arguments in
  \code{\dots} is either a point pattern object (of class \code{"ppp"})
  or a structure containing (at least) the elements
  \code{x} and \code{y}. The point patterns are not required to have
  the same window of observation.

  The window for the superimposed point pattern is
  specified by \code{W}. Its default value is
  the union of the windows of all the point patterns.

  If any of the arguments is a \emph{marked} point pattern,
  then all of them must be.
  In that case, the result is also a marked point pattern.

  If the arguments are given in the form \code{name=value},
  then the names will be used
  as marks attached to the corresponding points. See the last line
  in the Examples.
}
\seealso{
  \code{\link{concatxy}},
  \code{\link{quadscheme}}.

  See \code{\link{superimposePSP}} for superimposing patterns of line segments.
}
\examples{
  dat <- runifrect(30)
  xy <- list(x=runif(10),y=runif(10))
  new <- superimpose(dat, xy)
  \dontrun{plot(new)}

  # how to make a 2-type point pattern with types "a" and "b"
  u <- superimpose(a = rpoispp(10), b = rpoispp(20))

  # how to make a 2-type point pattern with types 1 and 2
  u <- superimpose("1" = rpoispp(10), "2" = rpoispp(20))
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{manip}
