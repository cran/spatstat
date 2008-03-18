\name{model.matrix.ppm}
\alias{model.matrix.ppm}
\title{Extract Design Matrix from Point Process Model}
\description{
  Given a point process model that has been
  fitted to spatial point pattern data,
  this function extracts the design matrix of the model.
}
\usage{
   \method{model.matrix}{ppm}(object, ...)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"}.)
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
}
\details{
  This command is a method for the generic function
  \code{\link{model.matrix}}. It extracts the design matrix of a
  spatial point process model (class \code{"ppm"}).

  More precisely, this command extracts
  the design matrix of the generalised linear model associated with
  a spatial point process model.
  
  The \code{object} must be a fitted point process model
  (object of class \code{"ppm"}) fitted to spatial point pattern data.
  Such objects are produced by the model-fitting
  function \code{\link{ppm}}.

  This function \code{model.matrix.ppm}
  extracts the model matrix for the GLM.

  The result is a matrix, with one row for every quadrature point
  in the fitting procedure, and one column for every constructed
  covariate in the design matrix. 

  The quadrature points themselves can be extracted using
  \code{\link{quad.ppm}}.
}
\value{
  A matrix. Rows of the matrix correspond to quadrature points
  in the fitting procedure. Columns are covariates in the model.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{model.images}},
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{quad.ppm}},
  \code{\link{residuals.ppm}},
  \code{\link{model.matrix}}
}

\examples{
   data(cells)
   fit <- ppm(cells, ~x)
   model.matrix(fit)
   # matrix with two columns: '(Intercept)' and 'x'
}
\keyword{spatial}
\keyword{models}
