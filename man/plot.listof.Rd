\name{plot.listof}
\alias{plot.listof}
\title{Plot a List of Things}
\description{
  Plots a list of things
}
\usage{
  \method{plot}{listof}(x, \dots, main, arrange=TRUE)
}
\arguments{
  \item{x}{
    An object of the class \code{"listof"}.
    Essentially a list of objects.
  }
  \item{\dots}{
    Arguments passed to \code{\link{plot}} when generating each
    plot panel.
  }
  \item{main}{
    Overall heading for the plot.
  }
  \item{arrange}{
    Logical flag indicating whether to plot the objects
    side-by-side. Set this to \code{FALSE} if you want to control
    the arrangement of the plots yourself.
  }
}
\value{
  Null.
}
\details{
  This is the \code{plot} method for the class \code{"listof"}.

  An object of class \code{"listof"} (defined in the base R package) represents
  a list of objects, all belonging to a common class.
  The base R package defines a method for printing these objects,
  \code{\link{print.listof}}, but does not define a method for \code{plot}.
  So here we have provided a method for \code{plot}.
  
  In the \pkg{spatstat} package, the function
  \code{\link{density.splitppp}} produces an object of
  class \code{"listof"}, essentially a list of pixel images.
  These images can be plotted in a nice arrangement
  using \code{plot.listof}. See the Example.
}
\seealso{
  \code{\link{print.listof}},
  \code{\link{density.splitppp}}
}
\examples{
# Multitype point pattern
 data(amacrine)
 plot(density(split(amacrine)))
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{hplot}
