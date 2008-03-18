\name{dummy.ppm}
\alias{dummy.ppm}
\title{Extract Dummy Points Used to Fit a Point Process Model}
\description{
  Given a fitted point process model,
  this function extracts the `dummy points' of the
  quadrature scheme used to fit the model.
}
\usage{
  dummy.ppm(object)
}
\arguments{
  \item{object}{
    fitted point process model (an object of class \code{"ppm"}).
  }
}
\value{
  A point pattern (object of class \code{"ppp"}).
}
\details{
  An object of class \code{"ppm"} represents a point process model
  that has been fitted to data. It is typically produced by
  the model-fitting algorithm \code{\link{ppm}}.
  
  The \code{\link{ppm}} algorithm approximates the pseudolikelihood
  integral by a sum over a finite set of quadrature points,
  which is constructed by augmenting
  the original data point pattern by a set of ``dummy'' points.
  The fitted model object returned by \code{\link{ppm}}
  contains complete information about this quadrature scheme.
  See \code{\link{ppm}} or \code{\link{ppm.object}} for further
  information.
  
  This function extracts the dummy points of the quadrature scheme.
  A typical use of this function would be to count the number of dummy
  points, to gauge the accuracy of the approximation to the
  exact pseudolikelihood. 

  See \code{\link{ppm.object}} for a list of all operations that can be
  performed on objects of class \code{"ppm"}.
}
\seealso{
  \code{\link{ppm.object}},
  \code{\link{ppp.object}},
  \code{\link{ppm}}
}
\examples{
 data(cells)
 fit <- ppm(cells, ~1, Strauss(r=0.1))
 X <- dummy.ppm(fit)
 X$n
 # this is the number of dummy points in the quadrature scheme
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{utilities}
\keyword{models}
