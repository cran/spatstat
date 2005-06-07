\name{rtoro}
\alias{rtoro}
\title{Random Toroidal Shift of Point Pattern}
\description{
  Applies a random shift simultaneously to all the points of a point pattern,
  or to selected sub-patterns, with wraparound at the
  borders of the window.
}
\usage{
   rtoro(X, which=NULL)
}
\arguments{
  \item{X}{Point pattern to be subjected to random shift.
    An object of class \code{"ppp"} or \code{"splitppp"}.
  }
  \item{which}{
    Identifies which components of the pattern will be shifted.
    Vector of marks (if \code{X} is a multitype point pattern)
    or subset of \code{names(X)} (if \code{X} is a \code{"splitppp"}
    object). Each component will be shifted independently.
  }
}
\value{
  An object of the same type as \code{X}.
}
\details{
  This operation randomises the locations of the points
  in a point pattern. It is useful for some hypothesis testing
  purposes.

  A `shift' is a simultaneous, parallel displacement of all points
  by the same displacement vector
  \eqn{(\Delta x, \Delta y)}{(delta.x, delta.y)}.
  Shifting a point pattern could leave some points outside the boundary of
  the window of observation, so we may (if the window is a rectangle)
  treat opposite edges of the window as identical, so that a point
  which disappears off the right-hand edge will re-appear at the
  left-hand edge. This is called a ``toroidal shift'' because it makes the
  rectangle topologically equivalent to the surface of a torus
  (doughnut). 

  The argument \code{X} may be either
  \itemize{
    \item
    a point pattern
    (an object of class \code{"ppp"})
    \item
    an object of class \code{"splitppp"}
    (basically a list of point patterns, obtained from
    \code{\link{split.ppp}}).
  }
  
  If \code{which=NULL}, then
  \itemize{
    \item
    If \code{X} is a point pattern object,
    all the points of the pattern
    are shifted simultaneously by the same displacement vector.
    \item
    If \code{X} is an object of class \code{"splitppp"}
    then each of the component point patterns is randomly shifted by a
    different displacement vector.
  }
  
  The argument \code{which} identifies which sub-patterns of
  \code{X} will be subjected to random shifts. Different sub-patterns
  will be shifted by different displacement vectors.
  If \code{which} is not \code{NULL}, then
  \itemize{
    \item
    If \code{X} is a point pattern object,
    it should be a multitype pattern (i.e. \code{X$marks} must be a
    factor). Then \code{which} should be a vector of types
    (i.e. a vector containing one or more levels of \code{X$marks}).
    For each value \code{a} in \code{which}, the sub-pattern of
    \code{X} consisting of points of type \code{a} only will be
    extracted, and subjected to a random toroidal shift.
    A different shift will be applied to different sub-patterns.
    The points not selected by \code{which} will be unchanged.
    \item
    If \code{X} is an object of class \code{"splitppp"},
    then \code{which} can be any valid subset index for the list \code{X}.
    For example \code{which} could be a vector of names of components
    in the list \code{X}. Each selected component point pattern will be shifted
    by a different displacement vector.
  }
  
  In this implementation, every random displacement vector is
  generated with a uniform distribution over the range of 
  possible displacement vectors. That is, the coordinates
  \eqn{\Delta x, \Delta y}{delta.x, delta.y} are independent
  and uniformly distributed.

}
\examples{
   data(amacrine)

   # point patterns:
   # shift all points simultaneously
   X <- rtoro(amacrine)
   # shift "on" and "off" points separately
   X <- rtoro(amacrine, which=amacrine$marks)
   # shift "on" points and leave "off" points fixed
   X <- rtoro(amacrine, which="on")

   # splitppp objects:
   Y <- split(amacrine)
   # shift "on" and "off" points separately
   Z <- rtoro(Y)
   # shift "on" points and leave "off" points fixed
   Z <- rtoro(Y, "on")
   # shift all points simultaneously
   Z <- split(rtoro(superimpose(Y)))
   
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
