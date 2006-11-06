\name{print.ppm}
\alias{print.ppm}
\title{Print a Fitted Point Process Model}
\description{
  Default \code{print} method for a fitted point process model.
}
\usage{
 \method{print}{ppm}(x,\dots)
}
\arguments{
  \item{x}{
    A fitted point process model, typically obtained from
    the model-fittingg algorithm \code{\link{ppm}}.
    An object of class \code{"ppm"}.
  }
  \item{\dots}{Ignored.}
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
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{plot.ppm}},
  \code{\link{predict.ppm}}
}
\examples{
 \dontrun{
 data(cells)
 Q <- quadscheme(cells)
 m <- ppm(Q, ~1, Strauss(0.05))
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
\keyword{print}
\keyword{models}

