\name{quad.ppm}
\alias{quad.ppm}
\title{Extract Quadrature Scheme Used to Fit a Point Process Model}
\description{
  Given a fitted point process model,
  this function extracts the 
  quadrature scheme used to fit the model.
}
\usage{
  quad.ppm(object)
}
\arguments{
  \item{object}{
    fitted point process model (an object of class \code{"ppm"}).
  }
}
\value{
  A quadrature scheme (object of class \code{"quad"}).
}
\details{
  An object of class \code{"ppm"} represents a point process model
  that has been fitted to data. It is typically produced by
  the maximum pseudolikelihood fitting algorithm \code{\link{mpl}}.
  
  The \code{\link{mpl}} algorithm approximates the pseudolikelihood
  integral by a sum over a finite set of quadrature points,
  which is constructed by augmenting
  the original data point pattern by a set of ``dummy'' points.
  The fitted model object returned by \code{\link{mpl}}
  contains complete information about this quadrature scheme.
  See \code{\link{mpl}} or \code{\link{ppm.object}} for further
  information.
  
  This function extracts the quadrature scheme.
  A typical use of this function would be to inspect the quadrature scheme
  (points and weights) to gauge the accuracy of the approximation to the
  exact pseudolikelihood. 

  See \code{\link{ppm.object}} for a list of all operations that can be
  performed on objects of class \code{"ppm"}.
  See \code{\link{quad.object}} for a list of all operations that can be
  performed on objects of class \code{"quad"}.
}
\seealso{
  \code{\link{ppm.object}},
  \code{\link{quad.object}},
  \code{\link{mpl}}
}
\examples{
 require(spatstat)
 data(cells)
 fit <- mpl(cells, ~1, Strauss(r=0.1), rbord=0.1)
 Q <- quad.ppm(fit)
 \dontrun{plot(Q)}
 Q$data$n
 Q$dummy$n
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
