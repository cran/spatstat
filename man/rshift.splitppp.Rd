\name{rshift.splitppp}
\alias{rshift.splitppp}
\title{Randomly Shift a List of Point Patterns}
\description{
  Randomly shifts each point pattern in a list of point patterns.
}
\synopsis{
   \method{rshift}{splitppp}(X, \dots, which=seq(X))
}
\usage{
   \method{rshift}{splitppp}(X, \dots, radius, width, height, edge="torus", clip,
                        which)
                        
}
\arguments{
  \item{X}{
    An object of class \code{"splitppp"}.
    Basically a list of point patterns.
  }
  \item{\dots}{
    Ignored.
  }
  \item{radius,width,height}{
    Optional parameters controlling the generation of the
    random shift vector. See Details.
  }
  \item{edge}{
    String indicating how to deal with edges of the patterns.
    Options are \code{"torus"}, \code{"erode"} and \code{"none"}.
    See Details.
  }
  \item{which}{
    Optional. Identifies which patterns will be shifted,
    while other patterns are not shifted.
    Any valid subset index for \code{X}.
  }
  \item{clip}{
    Optional. Window to which the final point patterns should be clipped.
  }
}
\value{
  Another object of class \code{"splitppp"}.
}
\details{
  This operation applies a random shift to each of the
  point patterns in the list \code{X}.
  
  The function \code{\link{rshift}} is generic.
  This function \code{rshift.splitppp}
  is the method for objects of class \code{"splitppp"}, which are
  essentially lists of point patterns, created by the function
  \code{\link{split.ppp}}.

  By default, every pattern in the list \code{X} will be shifted.
  The argument \code{which} indicates that only some of the patterns
  should be shifted, while other groups should be left unchanged.
  \code{which} can be any valid subset index for \code{X}.

  Each point pattern in the list \code{X}
  (or each pattern in \code{X[which]})
  is shifted by a random displacement vector.
  The shifting is performed by \code{\link{rshift.ppp}}.

  See the help page for \code{\link{rshift.ppp}}
  for details of the other arguments.
}
\seealso{
  \code{\link{rshift}},
  \code{\link{rshift.ppp}}
}
\examples{
   data(amacrine)
   Y <- split(amacrine)

   # random toroidal shift
   # shift "on" and "off" points separately
   X <- rshift(Y)

   # shift "on" points and leave "off" points fixed
   X <- rshift(Y, which="on")

   # maximum displacement distance 0.1 units
   X <- rshift(Y, radius=0.1)

   # shift with erosion
   X <- rshift(Y, radius=0.1, edge="erode")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
