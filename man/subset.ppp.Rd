\name{subset.ppp}
\alias{subset.ppp}
\alias{[.ppp}
\alias{[<-.ppp}
\title{Extract or Replace Subset of Point Pattern}
\description{
  Extract or replace a subset of a point pattern.
  Extraction has the effect of thinning the 
  points and/or or trimming the window.
}
\synopsis{
  \method{[}{ppp}(x, subset, window, drop, \dots)
  \method{[}{ppp}(x, subset, window) <- value
  subset.ppp(x, subset, window, drop, \dots)
}
\usage{
  \method{[}{ppp}(x, subset, window)
  \method{[}{ppp}(x, subset)
  \method{[}{ppp}(x, , window)

  \method{[}{ppp}(x, subset, window) <- value
  \method{[}{ppp}(x, subset) <- value
  \method{[}{ppp}(x, , window) <- value
}
\arguments{
  \item{x}{
    A two-dimensional point pattern.
    An object of class \code{"ppp"}.
  }
  \item{subset}{
    logical vector indicating which points should be retained.
  }
  \item{window}{
    window (an object of class \code{"owin"})
    delineating a subset of the original observation window. 
  }
  \item{value}{
    Replacement value for the subset. A point pattern.
  }
  \item{ }{ }
} 
\value{
  A point pattern (of class \code{"ppp"}).
}
\details{
  These functions extract a designated subset of a point pattern,
  or replace the designated subset with another point pattern.

  The function \code{[.ppp} is a method for \code{\link{[}} for the
  class \code{"ppp"}. It extracts a designated subset of a point pattern,
  either by ``\emph{thinning}''
  (retaining/deleting some points of a point pattern)
  or ``\emph{trimming}'' (reducing the window of observation
  to a smaller subregion and retaining only
  those points which lie in the subregion) or both.
 
  The pattern will be ``thinned''
  if  \code{subset} is specified. The points designated by \code{subset}
  will be retained. Here \code{subset} can be a numeric vector
  of positive indices (identifying the points to be retained),
  a numeric vector of negative indices (identifying the points
  to be deleted) or a logical vector of length equal to the number
  of points in the point pattern \code{x}. In the latter case,
  the points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  The pattern will be ``trimmed''
  if \code{window} is specified. This should
  be an object of class \code{\link{owin}} specifying a window of observation
  to which the point pattern \code{x} will be
  trimmed. The points of \code{x} lying inside the new
  \code{window} will be retained.
 
  Both ``thinning'' and ``trimming'' can be performed together.

  The function \code{[<-.ppp} is a method for \code{\link{[<-}} for the
  class \code{"ppp"}. It replaces the designated
  subset with the point pattern \code{value}.
  The subset of \code{x} to be replaced is designated by
  the arguments \code{subset} and \code{window} as above.

  The replacement point pattern \code{value} must lie inside the
  window of the original pattern \code{x}.
  The ordering of points in \code{x} will be preserved
  if the replacement pattern \code{value} has the same number of points
  as the subset to be replaced.  Otherwise the ordering is
  unpredictable.

  Use the function \code{\link{unmark}} to remove marks from a
  marked point pattern.

  Use the function \code{\link{split.ppp}} to select those points
  in a marked point pattern which have a specified mark.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{owin.object}},
  \code{\link{unmark}},
  \code{\link{split.ppp}},
  \code{\link{cut.ppp}}
}
\section{Warnings}{
  The function does not check whether \code{window} is a subset of
  \code{x$window}. Nor does it check whether \code{value} lies
  inside \code{x$window}.
}
\examples{
 data(longleaf)
 # Longleaf pines data
 \dontrun{
 plot(longleaf)
 }
 \testonly{
	longleaf <- longleaf[seq(1,longleaf$n,by=10)]
 }
 # adult trees defined to have diameter at least 30 cm
 adult <- (longleaf$marks >= 30)
 longadult <- longleaf[adult]
 # equivalent to: longadult <- subset.ppp(longleaf, subset=adult)
 \dontrun{
 plot(longadult)
 }
 # note that the marks are still retained.
 # Use unmark(longadult) to remove the marks
 
 # New Zealand trees data
 data(nztrees)
 \dontrun{
 plot(nztrees)          # plot shows a line of trees at the far right
 abline(v=148, lty=2)   # cut along this line
 }
 nzw <- owin(c(0,148),c(0,95)) # the subwindow
 # trim dataset to this subwindow
 nzsub <- nztrees[,nzw]
 # equivalent to: nzsub <- subset.ppp(nztrees, window=nzw)
 \dontrun{
 plot(nzsub)
 }

 # Redwood data
 data(redwood) 
 \dontrun{
 plot(redwood)
 }
 # Random thinning: delete 60\% of data
 retain <- (runif(redwood$n) < 0.4)
 thinred <- redwood[retain]
 \dontrun{
 plot(thinred)
 }
 # Scramble 60\% of data
 modif <- (runif(redwood$n) < 0.6)
 scramble <- function(x) { runifpoint(x$n, x$window) }
 redwood[modif] <- scramble(redwood[modif])

 # Lansing woods data - multitype points
 data(lansing)
 \testonly{
    lansing <- lansing[seq(1, lansing$n, length=100)]
 }
 # hickory trees only
 hick <- lansing[lansing$marks == "hickory", ]
 # still a marked pattern -- remove marks
 hick <- unmark(hick)

 # scramble the hickories
 lansing[lansing$marks == "hickory"] <- scramble(hick)
  
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
