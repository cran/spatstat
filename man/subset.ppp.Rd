\name{subset.ppp}
\alias{subset.ppp}
\alias{"[.ppp"}
\title{Extract Subset of Point Pattern}
\description{
  Extract a subset of a point pattern, either by thinning the 
  points or by trimming the window.
}
\synopsis{
  subset.ppp(x, subset, window, drop, \dots)
}
\usage{
  subset.ppp(x, subset, window)
  x[subset]
  x[ , window]
  x[subset,window]
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
} 
\value{
  A point pattern (of class \code{"ppp"}).
}
\details{
  This function takes performs either ``\emph{thinning}''
  (retaining/deleting some points of a point pattern)
  or ``\emph{trimming}'' (reducing the window of observation
  to a smaller subregion and retaining only
  those points which lie in the subregion). 
 
  The pattern will be ``thinned''
  if  \code{subset} is specified. This should
  be a logical vector of length equal to the number of points in the 
  point pattern \code{x}.
  The points \code{(x$x[i], x$y[i])} for which 
  \code{subset[i]=TRUE} will be retained, and the others
  will be deleted.
 
  The pattern will be ``trimmed''
  if \code{window} is specified. This should
  be an object of class \code{\link{owin}} specifying a window of observation
  to which the point pattern \code{x} will be
  trimmed. The points of \code{x} lying inside the new
  \code{window} will be retained.
 
  Both ``thinning'' and ``trimming'' can be performed together.

  Use the function \code{\link{unmark}} to remove marks from a
  marked point pattern.
}
\seealso{
  \code{\link{ppp.object}},
  \code{\link{owin.object}},
  \code{\link{unmark}}
}
\section{Warnings}{
  The function does not check whether \code{window} is a subset of
  \code{x$window}.
}
\examples{
 \dontrun{library(spatstat)}

 data(longleaf)
 # Longleaf pines data
 \dontrun{
 plot(longleaf)
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

 # Lansing woods data - multitype points
 data(lansing)
 # hickory trees only
 hick <- lansing[lansing$marks == "hickory", ]
 # still a marked pattern -- remove marks
 hick <- unmark(hick)
 
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
