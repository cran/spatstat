\name{predict.ppm}
\alias{predict.ppm}
\title{Prediction from a Fitted Point Process Model}
\description{
Given a fitted point process model obtained by \code{\link{mpl}},	
evaluate the spatial trend and the conditional intensity of the model
at new locations.
}
\synopsis{
predict.ppm(object, newdata, nx = 40, ny = NULL, type="trend",\dots)
}
\usage{
predict.ppm(object, nx = 40, ny = NULL, type="trend")
predict.ppm(object, newdata, type="trend")
}
\arguments{
  \item{object}{
    A fitted point process model, typically obtained from
    the maximum pseudolikelihood algorithm \code{\link{mpl}}. An object of
    class \code{"ppm"} (see \code{\link{ppm.object}}).
  }
  \item{newdata}{
    optional data frame giving the values of any spatial covariates
    for the locations at which the trend and conditional intensity should be
    computed. (Incompatible with \code{nx,ny})
  }
  \item{nx}{predictions are to be made at 
    an (\code{nx} by \code{ny}) rectangular grid of positions in the 
    original window of observation.
    (Incompatible with \code{newdata})
  }
  \item{ny}{
    see \code{nx} above; \code{ny} defaults to the same value as \code{nx}
  }
  \item{type}{
    character string.
    Indicates which property of the fitted model should be predicted.
    Options are \code{"trend"} for the spatial trend, and 
    \code{"cif"} or \code{"lambda"} for the conditional intensity.
  }
}
\value{
  \emph{If \code{newdata} is given:}
  a vector of predicted values for the spatial locations
  (and marks and covariates) given in the rows of \code{newdata}.

  \emph{If \code{newdata} is not given:}
  a list with components
  \code{x}, \code{y}, \code{z} giving the
  coordinates and predicted values for a rectangular grid of points.
  If \code{object} is an unmarked point process,
  \code{z} is a matrix such that \code{z[i,j]} 
  gives the predicted value at location \code{(x[i],y[j])};
  for marked point processes,
  \code{z} is a three-dimensional array 
  such that \code{z[i,j,m]} gives the predicted value for a
  point at location \code{(x[i],y[j])} with mark \code{m}.
  
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
  maximum pseudolikelihood algorithm \code{\link{mpl}}. This
  returns an object of class \code{"ppm"} representing 
  the fitted point process model (see \code{\link{ppm.object}}).
  The parameter estimates in this fitted model can be read off 
  simply by printing the \code{ppm} object.
  The spatial trend and conditional intensity of the 
  fitted model are evaluated using this function \code{predict.ppm}.
 
  This is feasible in such generality because the Berman-Turner-Baddeley method
  (Berman and Turner, 1992; Baddeley and Turner, 2000) 
  reduces maximum pseudolikelihood estimation to the fitting
  of a GLM or GAM, and because Splus/\R contains \code{predict} methods for
  \code{glm} and \code{gam} objects. The current implementation
  invokes \code{\link{predict.glm}}.
 
  The default action is to create a rectangular 40 by 40 grid of points
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
 
  The argument \code{newdata}, if given, is a data frame 
  giving the values of any spatial covariates
  at the locations where the trend and conditional intensity should be
  computed. 
  If the trend formula in the fitted model 
  involves spatial covariates (other than
  the Cartesian coordinates \code{x}, \code{y})
  then \code{newdata} is required.
 
  If \code{newdata} is present then it must contain variables
  matching all the variable names featuring in the trend formula
  \bold{and} the Cartesian coordinates \code{x} and \code{y}
  and the mark values \code{marks}.
  [This is different from the role of the \code{data}
  argument in \code{\link{mpl}} which must not contain
  \code{x} and \code{y} variables.]

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
  \code{\link{mpl}},
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
  If you are passing the \code{newdata} argument
  and the function crashes,
  it is advisable to start by checking that all the conditions 
  listed above are satisfied.
}
\examples{
  \dontrun{library(spatstat)}
  data(cells)
  m <- mpl(cells, ~ polynom(x,y,2), Strauss(0.05), rbord=0.05)
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
 
 
