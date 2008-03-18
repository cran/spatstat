\name{model.images}
\alias{model.images}
\title{Compute Images of Constructed Covariates}
\description{
  For a point process model fitted to spatial point pattern data,
  this function computes pixel images of the covariates
  in the design matrix.
}
\usage{
  model.images(object, W = as.owin(object), ...)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"}.)
  }
  \item{W}{
    A window (object of class \code{"owin"}) in which the
    images should be computed. Defaults to the window
    in which the model was fitted.
  }
  \item{\dots}{
    Other arguments (such as \code{na.action}) passed to
    \code{\link{model.matrix.lm}}.
  }
}
\details{
  This command is similar to \code{\link{model.matrix.ppm}} except
  that it computes pixel images of the covariates,
  instead of computing the covariate values at certain points only.

  The \code{object} must be a fitted spatial point process model
  (object of class \code{"ppm"}) produced by the model-fitting
  function \code{\link{ppm}}.

  The spatial covariates required by the model-fitting procedure
  are computed at every pixel location in the window \code{W}.

  Note that the spatial covariates computed here
  are not the original covariates that were supplied when fitting the
  model. Rather, they are the covariates that actually appear in the
  loglinear representation of the (conditional) intensity
  and in the columns of the design matrix. For example, they might include
  dummy or indicator variables for different levels of a factor,
  depending on the contrasts that are in force.

  The pixel resolution is determined by \code{W} 
  if \code{W} is a mask (that is \code{W$type = "mask"}).
  Otherwise, the pixel resolution is determined by
  \code{\link{spatstat.options}}.

  The result is a named list of pixel images (objects of class
  \code{"im"}) containing the values of the spatial covariates.
  The names of the list elements are the names of the covariates
  determined by \code{\link{model.matrix.lm}}.
  
  The result is also of class \code{"listof"} so that it can
  be plotted immediately.
}
\value{
  An object of class \code{"listof"} consisting of 
  a named list of pixel images (objects of class \code{"im"}).
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\seealso{
  \code{\link{model.matrix.ppm}},
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{im}},
  \code{\link{im.object}},
  \code{\link{plot.listof}},
  \code{\link{spatstat.options}}
}

\examples{
   data(cells)
   fit <- ppm(cells, ~x)
   model.images(fit)
   fit2 <- ppm(cells, ~cut(x,3))
   model.images(fit2)
}
\keyword{spatial}
\keyword{models}
