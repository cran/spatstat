\name{coords}
\Rdversion{1.1}
\alias{coords}
\alias{coords.ppp}
\alias{coords.ppx}
\title{
  Extract Coordinates of a Spatial or Spatiotemporal Point Pattern
}
\description{
  Given any kind of spatial or space-time point pattern,
  this function extracts the
  (space and/or time) coordinates of the points
  and returns them as a data frame.
}
\usage{
  coords(x, ...)
  \method{coords}{ppp}(x, ...)
  \method{coords}{ppx}(x, ..., spatial = TRUE, temporal = TRUE)
}
\arguments{
  \item{x}{
    A point pattern: either a two-dimensional point pattern
    (object of class \code{"ppp"}), a three-dimensional point pattern
    (object of class \code{"pp3"}), or a 
    general multidimensional space-time point pattern 
    (object of class \code{"ppx"}).
  }
  \item{\dots}{
    Further arguments passed to methods.
  }
  \item{spatial,temporal}{
    Logical values indicating whether to extract spatial and temporal
    coordinates, respectively. The default is to return both
    spatial and temporal coordinates.
  }
}
\details{
  The function \code{coords} is generic, with methods for
  the classes \code{"ppp"}) and \code{"ppx"}.
  An object of class \code{"pp3"} also inherits from \code{"ppx"} and
  is handled by the method for \code{"ppx"}.
}
\value{
  A \code{data.frame} with one row for each point, containing the
  coordinates.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{ppx}},
  \code{\link{pp3}},
  \code{\link{ppp}},
  \code{as.hyperframe.ppx},
  \code{as.data.frame.ppx}.
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),t=runif(4))
   X <- ppx(data=df, temporal="t")
   coords(X)
   coords(X, temporal=FALSE)
}
\keyword{spatial}
\keyword{manip}
