\name{plot.ppp}
\alias{plot.ppp}
\title{plot a Spatial Point Pattern}
\description{
  Plot a two-dimensional spatial point pattern
}
\usage{
 plot.ppp(x, main, \dots, chars, use.marks=TRUE, add=FALSE)
}
\arguments{
  \item{x}{
    The spatial point pattern to be plotted.
    An object of class \code{"ppp"},
    or data which can be converted into 
    this format by \code{\link{as.ppp}()}.
  }
  \item{main}{
    text to be displayed as a title above the plot.
  }
  \item{\dots}{
    extra arguments passed to the generic \code{\link{plot}} function.
  }
  \item{chars}{
    the plotting characters which should be used to represent
    the points with different marks. See below.
  }
  \item{use.marks}{
    logical flag; if \code{TRUE}, plot points using a different
    plotting symbol for each mark;
    if \code{FALSE}, only the locations of the points will be plotted,
    using \code{\link{points}()}.
  }
  \item{add}{
    logical flag; if \code{TRUE},
    just the points are plotted, over the existing plot.
    A new plot is not created, and
    the window is not plotted.
  }
}
\value{
  \code{NULL}, or a vector giving the correspondence between 
  mark values and plotting characters.
}
\details{
  This is the \code{plot} method for the class \code{"ppp"}
  (see \code{\link{ppp.object}} for details of this class).
  
  It invokes \code{\link{plot.owin}} to display the observation window
  \code{x$window}.
  Then the locations of the points are
  plotted using \code{\link{points}()}. 
 
  If the points do not have marks (i.e. \code{x$marks} is absent)
  or if \code{use.marks = FALSE},
  then the locations of the points will be plotted 
  using a single plot character. This character can be changed
  by the extra argument \code{pch};  see the examples.
 
  If the points have marks, then each mark value will be 
  represented by a different plot character. This action can be overridden
  by setting \code{use.marks = FALSE}.
 
  For marked point patterns the default 
  is to represent the \eqn{i}th unique mark value by
  \code{points(..., pch=i)}. 
  If the vector \code{chars} is supplied, it should have the same length
  as \code{sort(unique(x$marks))}, and then the \eqn{i}th smallest mark
  will be plotted using \code{points(..., pch=chars[i])}. 
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{plot}},
  \code{\link{points}},
  \code{\link{plot.owin}}
}
\examples{
   \dontrun{
   library(spatstat)

   data(cells)
   plot(cells)

   # multitype 
   data(lansing)
   plot(lansing)

   # marked by a real number
   data(longleaf)
   plot(longleaf)

   # just plot the points
   plot(longleaf, use.marks=F)
   plot(unmark(longleaf)) # equivalent

   # use selected plot characters
   lma <- sort(unique(longleaf$marks))
   adult <- ifelse(lma >= 30, "o", "+")
   plot(longleaf, chars = adult)

  }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
 
