\name{superimpose}
\alias{superimpose}
\title{Superimpose Several Point Patterns}
\description{
  Superimpose any number of point patterns.
}
\usage{
  superimpose(\dots)
}
\arguments{
  \item{\dots}{
    Any number of arguments, each of which represents a point pattern.
    Each argument must be either a point pattern
    object (of class \code{"ppp"}) or a structure
    containing elements \code{x} and \code{y}.
  }
}
\value{
  A point pattern (object of class \code{"ppp"}) representing the
  superposition (union) of all the point patterns.
}
\details{
  This function is used to superimpose two or more point patterns
  (see also \code{\link{concatxy}}). 
  
  It assumes that each of the arguments in
  \code{\dots} is either a point pattern object (of class \code{"ppp"})
  or a structure containing (at least) the elements
  \code{x} and \code{y}. All the point patterns are assumed to have
  the same window of observation, and (if they are marked point
  patterns) the same space of possible marks.

  The first argument must be a point pattern
  object, and it determines the observation window and mark space
  for the result.

  To superimpose point patterns consisting only of \code{x,y} coordinates,
  use \code{\link{concatxy}}.
}
\seealso{
  \code{\link{concatxy}},
  \code{\link{quadscheme}}
}
\examples{
  dat <- runifrect(30)
  xy <- list(x=runif(10),y=runif(10))
  new <- superimpose(dat, xy)
  \dontrun{plot(new)}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
