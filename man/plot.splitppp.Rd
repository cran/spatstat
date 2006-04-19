\name{plot.splitppp}
\alias{plot.splitppp}
\title{Plot a List of Point Patterns}
\description{
  Plots a list of point patterns.
}
\usage{
  \method{plot}{splitppp}(x, \dots, arrange=TRUE)
}
\arguments{
  \item{x}{
    A named list of point patterns,
    typically obtained from \code{\link{split.ppp}}.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot.ppp}} which control the
    appearance of each plot panel.
  }
  \item{arrange}{
    Logical flag indicating whether to plot the point patterns
    side-by-side. Set this to \code{FALSE} if you want to control
    the arrangement of the plots yourself.
  }
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"splitppp"}.
  It is typically used to plot the result of the function
  \code{\link{split.ppp}} but it may also be used to plot any list
  of point patterns created by the user.

  The argument \code{x} should be a named list of point patterns
  (objects of class \code{"ppp"}, see \code{\link{ppp.object}}).
  Each of these point patterns will be plotted in turn
  using \code{\link{plot.ppp}}. 
}
\seealso{
  \code{\link{split.ppp}},
  \code{\link{plot.ppp}},
  \code{\link{ppp.object}}
}
\examples{
# Multitype point pattern
 data(amacrine)
 plot(split(amacrine))
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
