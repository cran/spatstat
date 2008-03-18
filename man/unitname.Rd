\name{unitname}
\alias{unitname}
\alias{unitname.owin}
\alias{unitname.psp}
\alias{unitname.ppp}
\alias{unitname.im}
\alias{unitname<-}
\alias{unitname<-.ppp}
\alias{unitname<-.psp}
\alias{unitname<-.owin}
\alias{unitname<-.im}
\title{Name for Unit of Length}
\description{
  Inspect or change the name of the unit of length
  in a spatial dataset.
}
\usage{
unitname(x)
\method{unitname}{ppp}(x)
\method{unitname}{psp}(x)
\method{unitname}{owin}(x)
\method{unitname}{im}(x)
unitname(x) <- value
\method{unitname}{ppp}(x) <- value
\method{unitname}{psp}(x) <- value
\method{unitname}{owin}(x) <- value
\method{unitname}{im}(x) <- value
}
\arguments{
  \item{x}{A spatial dataset.
    Either a point pattern (object of class \code{"ppp"}),
    a line segment pattern (object of class \code{"psp"}),
    a window (object of class \code{"owin"}) or
    a pixel image (object of class \code{"im"}).
  }
  \item{value}{
    Name of the unit of length. See Details.
  }
}
\details{
  Spatial datasets in the \pkg{spatstat} package
  may include the name of the unit of length. This name is used
  when printing or plotting the dataset, and in some other
  applications. 

  \code{unitname(x)} extracts this name,
  and \code{unitname(x) <- value} sets the name to \code{value}.

  A valid name is either
  \itemize{
    \item a single character string
    \item a vector of two character strings giving the
    singular and plural forms of the unit name
    \item a list of length 3, containing two character strings
    giving the singular and plural forms of the basic unit,
    and a number specifying the multiple of this unit.
  }
}
\value{
  The return value of \code{unitname} is an object of class \code{"units"}
  containing the name of the unit of length in \code{x}. There are
  methods for \code{print} and \code{summary}.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{owin}},
  \code{\link{ppp}}
}

\examples{
  X <- runifpoint(20)

  # if the unit of length is 1 metre:
  unitname(X) <- c("metre", "metres")

  # if the unit of length is 6 inches:
  unitname(X) <- list("inch", "inches", 6)
}
\keyword{spatial}
\keyword{manip}

