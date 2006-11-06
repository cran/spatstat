\name{subset.fasp}
\alias{subset.fasp}
\alias{[.fasp}
\title{Extract Subset of Function Array}
\description{
  Extract a subset of a function array (an object of class
  \code{"fasp"}).
}
\usage{
  \method{[}{fasp}(x, I, J, drop,\dots)
}
\arguments{
  \item{x}{
    A function array. An object of class \code{"fasp"}.
  }
  \item{I}{
    any valid expression for a subset of the row indices of the array.
  }
  \item{J}{
    any valid expression for a subset of the column indices of the array.
  }
  \item{drop, \dots}{Ignored.}
} 
\value{
  A function array (of class \code{"fasp"}).
}
\details{
  A function array can be regarded as a matrix whose entries
  are functions. See \code{\link{fasp.object}} for an explanation of
  function arrays.

  This routine extracts a sub-array according to the usual
  conventions for matrix indexing.
}
\seealso{
  \code{\link{fasp.object}}
}
\examples{
 # Lansing woods data - multitype points with 6 types
 data(lansing)
 \testonly{
 # smaller dataset
   lansing <- lansing[ seq(1,lansing$n,by=45), ]
 }
 # compute 6 x 6 array of all cross-type K functions
 a <- alltypes(lansing, "K")

 # extract first three marks only
 b <- a[1:3,1:3]
 \dontrun{plot(b)}
 # subset pertaining to hickories
 h <- a[levels(lansing$marks) == "hickory", ]
 \dontrun{plot(h)}
}

\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
\keyword{manip}
