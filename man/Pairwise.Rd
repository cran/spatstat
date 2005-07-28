\name{Pairwise}
\alias{Pairwise}
\title{Generic Pairwise Interaction model}
\description{
Creates an instance of a pairwise interaction point process model
which can then be fitted to point pattern data.
}
\synopsis{
  Pairwise(pot, name, par, parnames, printfun)
}
\usage{
  Pairwise(pot, name)
  Pairwise(pot, name, par, parnames)
}
\arguments{
  \item{pot}{An S language function giving the user-supplied
    pairwise interaction potential.}
  \item{name}{Character string.}
  \item{par}{List of numerical values for irregular parameters}
  \item{parnames}{Vector of names of irregular parameters}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  Advanced use only.
  
  This code constructs a member of the
  pairwise interaction family \code{\link{pairwise.family}}
  with arbitrary pairwise interaction potential given by
  the user.

  The function \code{pot} must take as its first argument
  a matrix of interpoint distances, and evaluate the
  potential for each of these distances. The result must be
  either a matrix with the same dimensions as its input,
  or an array with its first two dimensions the same as its input
  (the latter case corresponds to a vector-valued potential).

  If irregular parameters are present, then the second argument
  to \code{pot} should be a vector of the same type as \code{par}
  giving those parameter values.  
}
\seealso{
  \code{\link{ppm}},
  \code{\link{pairwise.family}},
  \code{\link{ppm.object}}
}
\examples{
   #This is the same as StraussHard(r=0.7,h=0.2)
   strpot <- function(d,par) {
         r <- par$r
         h <- par$h
         value <- (d <= r)
         value[d < h] <- -Inf
         value
   }
   mySH <- Pairwise(strpot, "StraussHard", list(r=0.7,h=0.2),
           c("interaction distance r", "hard core distance c"))
   data(cells)
   ppm(cells, ~ 1, mySH)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}

