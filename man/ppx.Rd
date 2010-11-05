\name{ppx}
\Rdversion{1.1}
\alias{ppx}
\title{
  Multidimensional Space-Time Point Pattern
}
\description{
  Creates a multidimensional space-time point pattern
  with any kind of coordinates and marks.
}
\usage{
ppx(data, domain=NULL, spatial = NULL, temporal = NULL)
}
\arguments{
  \item{data}{
    The coordinates and marks of the points.
    A \code{data.frame} or \code{hyperframe}.
  }
  \item{domain}{
    Optional. The space-time domain containing the points.
    An object in some appropriate format, or \code{NULL}.
  }
  \item{spatial}{
    Character vector giving the names of the columns of
    \code{data} that should be interpreted as spatial coordinates.
  }
  \item{temporal}{
    Character vector giving the names of the columns of
    \code{data} that should be interpreted as temporal coordinates.
  }
}
\details{
  An object of class \code{"ppx"}
  represents a marked point pattern
  in multidimensional space and/or time. There may be any
  number of spatial coordinates, any number of temporal coordinates,
  and any number of mark variables. The individual marks may be
  atomic (numeric values, factor values, etc) or objects of any kind.

  The argument \code{data} should contain the coordinates and marks of
  the points. It should be a \code{data.frame} or more generally a
  \code{hyperframe} (see \code{\link{hyperframe}}) with
  one row of data for each point.

  Each column of \code{data} is either
  a spatial coordinate, a temporal coordinate, or a mark variable.
  The argument \code{spatial} determines which columns are interpreted
  as spatial coordinates: if it is present, \code{spatial} should be a
  character vector matching the names of columns of \code{data}. 
  Similarly the argument \code{temporal} determines which columns
  are interpreted as temporal coordinates.
  
  By default,  columns of numerical data are assumed to represent
  spatial coordinates.
}
\value{
  An object of class \code{"ppx"}.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{pp3}},
  \code{\link{print.ppx}} 
}
\examples{
   df <- data.frame(x=runif(4),y=runif(4),z=runif(4))
   X <- ppx(data=df, t="z")
   X

   val <- runif(4)
   E <- lapply(val, function(s) { rpoispp(s) })
   hf <- hyperframe(t=val, e=as.listof(E))
   Z <- ppx(data=hf, domain=c(0,1))
   Z
}
\keyword{spatial}
\keyword{datagen}
