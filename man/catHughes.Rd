\name{catHughes}
\alias{catHughes}
\title{Hughes' Cat Retina Data}
\description{
Austin Hughes' data: a point pattern 
of ganglion cells in the retina of a cat.
A marked point pattern.
} 
\format{
  An object of class \code{"ppp"}
  representing the point pattern of cell locations.
  Entries include
  \tabular{ll}{
    \code{x} \tab Cartesian \eqn{x}-coordinate of cell \cr
    \code{y} \tab Cartesian \eqn{y}-coordinate of cell \cr
    \code{marks} \tab factor with levels 0 and 1 indicating
      ``off'' and ``on'' cells respectively.
    }
  See \code{\link{ppp.object}} for details of the format.
}
\usage{data(catHughes)}
\source{Peter Diggle, personal communication}
\section{Notes}{
Austin Hughes' data: a point pattern 
of ganglion cells in the retina of a cat.
152 ``on'' cells and 142 ``off'' cells in a rectangle
with dimensions 1.6065 by 1.0 where the unit of measurement
is approximately 1 mm. 
}
\keyword{datasets}
