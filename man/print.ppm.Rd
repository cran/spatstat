\name{print.ppm}
\alias{print.ppm}
\title{Print a Fitted Point Process Model}
\description{
  Default \code{print} method for a fitted point process model.
}
\synopsis{
 print.ppm(x,\dots)
}
\usage{
 print.ppm(x)
 x
}
\arguments{
  \item{x}{
    A fitted point process model, typically obtained from
    the maximum pseudolikelihood algorithm \code{\link{mpl}}.
    An object of class \code{"ppm"}.
  }
}
\value{
  none.
}
\details{
  This is the \code{print} method for the class \code{"ppm"}.
  It prints information about the fitted model in a sensible format.

  See \code{\link{ppm.object}} for details of the class \code{"ppm"}.
}
\seealso{
  \code{\link{mpl}},
  \code{\link{ppm.object}},
  \code{\link{plot.ppm}},
  \code{\link{predict.ppm}}
}
\examples{
 \dontrun{
 library(spatstat)
 data(cells)
 Q <- quadscheme(cells)
 m <- mpl(Q, ~1, Strauss(0.05))
 m
 }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
