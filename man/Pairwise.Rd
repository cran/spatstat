\name{Pairwise}
\alias{Pairwise}
\title{Generic Pairwise Interaction model}
\description{
Creates an instance of a pairwise interaction point process model
which can then be fitted to point pattern data.
}
\usage{
  Pairwise(pot, name)
}
\arguments{
  \item{pot}{An S language function giving the user-supplied
    pairwise interaction potential.}
  \item{name}{Character string.}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  This is experimental. It constructs a member of the
  pairwise interaction family \code{\link{pairwise.family}}.
}
\seealso{
  \code{\link{mpl}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}

