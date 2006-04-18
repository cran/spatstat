\name{update.ppm}
\alias{update.ppm}
\title{Update a Fitted Point Process Model}
\description{
  \code{update} method for class \code{"ppm"}.
}
\usage{
  \method{update}{ppm}(object, \dots, fixdummy=TRUE)

}
\arguments{
  \item{object}{
    An existing fitted point process model,
    typically produced by \code{\link{ppm}}.
  }
  \item{\dots}{
    Arguments to be updated in the new call to \code{\link{ppm}}.
  }
  \item{fixdummy}{
    Logical flag indicating whether the quadrature scheme for the
    call to \code{\link{ppm}} should use the same set of dummy points
    as that in the original call.
  }
}
\details{
  This is a method for the generic function \code{\link{update}}
  for the class \code{"ppm"}. An object of class \code{"ppm"}
  describes a fitted point process model. See \code{\link{ppm.object}})
  for details of this class.

  \code{update.ppm} will modify the point process model
  specified by \code{object} according to the new arguments given,
  then re-fit it.
  The actual re-fitting is performed by the model-fitting
  function \code{\link{ppm}}.

  If you are comparing several model fits to the same data,
  or fits of the same model to different data, it is
  strongly advisable to use \code{update.ppm}
  rather than trying to fit them by hand.
  This is because \code{update.ppm} re-fits the model
  in a way which is comparable to the original fit.

  The arguments \code{...} are matched to the formal arguments
  of \code{\link{ppm}} as follows.

  First, all the \emph{named} arguments in \code{...} are matched
  with the formal arguments of \code{\link{ppm}}.
  Use \code{name=NULL} to remove the argument \code{name} from the
  call.

  Second, any \emph{unnamed} arguments in \code{...} are
  matched with formal arguments of \code{\link{ppm}} if the matching
  is obvious from the class of the object. Thus \code{...} may contain
  \itemize{
    \item
    exactly one argument of class \code{"ppp"} or \code{"quad"},
    which will be interpreted as the named argument \code{Q};
    \item
    exactly one argument of class \code{"formula"}, which will be
    interpreted as the named argument \code{trend};
    \item
    exactly one argument of class \code{"interact"}, which will be
    interpreted as the named argument \code{interaction};
    \item
    exactly one argument of class \code{"data.frame"}, which will be
    interpreted as the named argument \code{covariates}.
  }

  When \code{fixdummy=FALSE}, calling \code{update.ppm}
  is exactly the same as calling \code{ppm} with the updated
  arguments. However, the original and updated models
  are not strictly comparable (for example, their pseudolikelihoods
  are not strictly comparable) unless they used the same set of dummy
  points for the quadrature scheme. Setting \code{fixdummy=TRUE}
  ensures that the re-fitting will be performed using the same set
  of dummy points. This is highly recommended.
}
\value{
  Another fitted point process model (object of class \code{"ppm"}).
}
\examples{
  data(nztrees)
  data(cells)

  # fit the stationary Poisson process
  fit <- ppm(nztrees, ~ 1)

  # fit a nonstationary Poisson process
  fitP <- update(fit, trend=~x)
  fitP <- update(fit, ~x)

  # fit a stationary Strauss process
  fitS <- update(fit, interaction=Strauss(13))
  fitS <- update(fit, Strauss(13))

  # oops, forgot the edge correction
  fitS <- update(fitS, rbord=13)

  # re-fit the model to a subset
  # of the original point pattern
  nzw <- owin(c(0,148),c(0,95))
  nzsub <- nztrees[,nzw]
  fut <- update(fitS, Q=nzsub)
  fut <- update(fitS, nzsub)

  # WARNING: the point pattern argument is called 'Q'
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
