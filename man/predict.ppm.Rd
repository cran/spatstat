\name{predict.ppm}
\alias{predict.ppm}
\title{Prediction from a Fitted Point Process Model}
\description{
Given a fitted point process model obtained by \code{\link{ppm}},	
evaluate the spatial trend and the conditional intensity of the model
at new locations.
}
\usage{
   predict.ppm(object, window, ngrid=NULL, locations=NULL,
   covariates=NULL, type="trend", \dots)
}
\arguments{
  \item{object}{
    A fitted point process model, typically obtained from
    the model-fitting algorithm \code{\link{ppm}}. An object of
    class \code{"ppm"} (see \code{\link{ppm.object}}).
  }
  \item{window}{
    window delimiting the locations where predictions
    should be computed. Defaults to the window of the
    original data used to fit the model \code{object}.
  }
  \item{ngrid}{dimensions (either (\code{ngrid[1]} by \code{ngrid[2]}
    or \code{ngrid} by \code{ngrid}) of a rectangular grid of locations
    inside \code{window} where predictions should be computed.
    (Incompatible with \code{locations})
  }
  \item{locations}{
    data giving the \eqn{x,y} coordinates (and marks, if required)
    at which predictions should be computed. 
    (Incompatible with \code{ngrid})
  }
  \item{covariates}{
    Values of external covariates required by the model.
    Either a data frame or a list of images.
    See Details.
  }
  \item{type}{
    character string.
    Indicates which property of the fitted model should be predicted.
    Options are \code{"trend"} for the spatial trend, and 
    \code{"cif"} or \code{"lambda"} for the conditional intensity.
  }
  \item{\dots}{
    Ignored.
  }
}
\value{
  \emph{If \code{locations} is given:}
  a vector of predicted values for the spatial locations
  (and marks, if required) given in \code{locations}.

  \emph{If \code{locations} is not given:}
  If \code{object} is an unmarked point process,
  the result is an image object (of class \code{"im"}, see
  \code{\link{im.object}}) containing the predictions. 
  If \code{object} is a multitype point process,
  the result is a list of images, containing the predictions
  for each type at the same grid of locations.

  The ``predicted values'' are either values of the spatial trend
  (if \code{type="trend"}) or values of the conditional intensity
  (if \code{type="cif"} or \code{type="lambda"}).
}
\details{
  This function computes the spatial trend
  and the conditional intensity of a fitted spatial point process model.
  See Baddeley and Turner (2000) for explanation and examples.
 
  Given a point pattern dataset, we may fit
  a point process model to the data using the 
  model-fitting algorithm \code{\link{ppm}}. This
  returns an object of class \code{"ppm"} representing 
  the fitted point process model (see \code{\link{ppm.object}}).
  The parameter estimates in this fitted model can be read off 
  simply by printing the \code{ppm} object.
  The spatial trend and conditional intensity of the 
  fitted model are evaluated using this function \code{predict.ppm}.
 
  This is feasible in such generality because the Berman-Turner-Baddeley method
  (Berman and Turner, 1992; Baddeley and Turner, 2000) 
  reduces maximum pseudolikelihood estimation to the fitting
  of a GLM or GAM, and because \R contains \code{predict} methods for
  \code{glm} and \code{gam} objects. The current implementation
  invokes \code{\link{predict.glm}}.

  The default action is to create a rectangular grid of points
  in the observation window of the data point pattern, and evaluate
  the spatial trend and conditional intensity at these locations.
  Note that by ``spatial trend'' we mean the
  (exponentiated) first order potential
  and not the intensity of the process. [For example if we fit the
  stationary Strauss process with parameters
  \eqn{\beta}{beta} and \eqn{\gamma}{gamma},
  then the spatial trend is constant and equal to \eqn{\beta}{beta}. ]
  The conditional intensity \eqn{\lambda(u, X)}{lambda(u,X)} of the
  fitted model is evaluated at each required spatial location \eqn{u},
  with respect to the data point pattern \eqn{X}.

  If the argument \code{locations} is absent, then
  predictions are performed at an \code{ngrid} by \code{ngrid}
  pixel grid of locations in the window \code{window}. 
  The result of prediction will be a pixel image or images.
  
  If \code{locations} is present, then predictions
  will be performed at the locations given by
  this dataset. The result of prediction will be a
  vector of values.

  For an \emph{unmarked} point process model,
  \code{locations} is expected to specify the \eqn{x,y} coordinates
  of the prediction locations, in any format acceptable
  to \code{\link{xy.coords}}. For a marked point process model,
  \code{locations} must specify the \eqn{x,y} coordinates
  \bold{and} the \emph{marks} of the prediction locations;
  it should be a data frame or list containing entries named
  \code{x}, \code{y} and \code{marks}.
  
  The argument \code{covariates} gives the values of any spatial covariates
  at the prediction locations.
  If the trend formula in the fitted model 
  involves spatial covariates (other than
  the Cartesian coordinates \code{x}, \code{y})
  then \code{covariates} is required.

  The format and use of \code{covariates} is analogous to that of the
  argument of the same name in \code{\link{ppm}}.
  It is either a data frame or a list of images.
  If it is a data frame, then each row of \code{covariates}
  is assumed to contain covariate data for the location
  specified in the corresponding row of \code{locations}
  (which must be present).

  If \code{covariates} is a list of images, then
  the names of the entries should correspond to
  the names of covariates in the model formula \code{trend}.
  Each entry in the list must be an image object (of class \code{"im"},
  see \code{\link{im.object}}).
  The software will look up
  the pixel values of each image at the quadrature points.

  Note that if you only want to use prediction in order to
  generate a plot of the predicted values,
  it may be easier to use \code{\link{plot.ppm}} which calls
  this function and plots the results.
  
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Berman, M. and Turner, T.R. 
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41} (1992) 31--38.
}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{plot.ppm}},
  \code{\link{print.ppm}}
}
\section{Warnings}{
  The current implementation invokes \code{\link{predict.glm}}
  so that \bold{prediction is wrong} if the trend formula in
  \code{object} involves terms in \code{ns()},
  \code{bs()} or \code{poly()}.
  This is a weakness of \code{\link{predict.glm}} itself!
  
  Error messages may be very opaque,
  as they tend to come from deep in the workings of 
  \code{\link{predict.glm}}.
  If you are passing the \code{covariates} argument
  and the function crashes,
  it is advisable to start by checking that all the conditions 
  listed above are satisfied.
}
\examples{
  data(cells)
  m <- ppm(cells, ~ polynom(x,y,2), Strauss(0.05), rbord=0.05)
  trend <- predict(m, type="trend")
  \dontrun{
  image(trend)
  points(cells)
  }
  cif <- predict(m, type="cif")
  \dontrun{
  persp(cif)
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
 
 
