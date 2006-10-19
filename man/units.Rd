\name{units}
\alias{units}
\alias{units.owin}
\alias{units.ppp}
\alias{units.im}
\alias{units<-}
\alias{units<-.ppp}
\alias{units<-.owin}
\alias{units<-.im}
\title{Name for Unit of Length}
\description{
  Inspect or change the name of the unit of length
  in a spatial dataset.
}
\usage{
units(x)
\method{units}{ppp}(x)
\method{units}{owin}(x)
\method{units}{im}(x)
units(x) <- value
\method{units}{ppp}(x) <- value
\method{units}{owin}(x) <- value
\method{units}{im}(x) <- value
}
\arguments{
  \item{x}{A spatial dataset.
    Either a point pattern (object of class \code{"ppp"}),
    a window (object of class \code{"owin"}) or
    a pixel image (object of class \code{"im"}).
  }
  \item{value}{Character string, or vector of two character strings,
    containing the name of the unit of length.
  }
}
\details{
  Spatial datasets in the \pkg{spatstat} package
  may include the name of the unit of length. This name is used
  when printing or plotting the dataset, and in some other
  applications. 

  \code{units(x)} extracts this name as a character string,
  and \code{units(x) <- value} sets the name to \code{value}.

  A valid name is either a single character string,
  or a vector of two character strings giving the
  singular and plural forms.
}
\value{
  The value of \code{units(x)} is a character string
  containing the name of the unit of length in \code{x}.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\seealso{
  \code{\link{owin}},
  \code{\link{ppp}}
}

\examples{
  data(bei)
  units(bei) <- c("metre", "metres")
}
\keyword{spatial}
\keyword{manip}

