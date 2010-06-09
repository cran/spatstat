\name{model.frame.ppm}
\alias{model.frame.ppm}
\title{
Extract the Environment of a Point Process Model
}
\description{
  Given a fitted point process model, this function
  returns a data frame containing all the variables needed to
  fit the model using the Berman-Turner device.
}
\usage{
 \method{model.frame}{ppm}(formula, ...)
}
\arguments{
  \item{formula}{
    A fitted point process model (object of class \code{"ppm"}).
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{model.frame.glm}}.
  }
}
\details{
  The function \code{\link{model.frame}} is generic.
  This function is a method for \code{\link{model.frame}}
  for fitted point process models (objects of class \code{"ppm"}).

  The first argument should be a fitted point process model;
  it has to be named \code{formula} for consistency with the generic
  function.

  The result is a data frame containing all the variables used in
  fitting the model. The data frame has one row for each quadrature point
  used in fitting the model. The quadrature scheme can be extracted using
  \code{\link{quad.ppm}}.  
}
\value{
  A \code{data.frame} containing all the variables used in the
  fitted model, plus additional variables specified in \code{\dots}.
  It has an additional attribute \code{"terms"} containing information
  about the model formula. For details see \code{\link{model.frame.glm}}.
}
\references{
  Baddeley, A. and Turner, R. (2000)
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42}, 283--322.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{model.frame}},
  \code{\link{model.matrix.ppm}}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~x)
  mf <- model.frame(fit)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

