\name{update.ppm}
\alias{update.ppm}
\title{Update a Fitted Point Process Model}
\description{
  \code{update} method for class \code{"ppm"}.
}
\usage{
  update.ppm(object, \dots,
            Q, trend, interaction, covariates,
            correction, rbord, use.gam)
}
\arguments{
  \item{object}{
    An existing fitted point process model,
    typically produced by \code{\link{ppm}}.
  }
  \item{\dots}{
    Unnamed arguments - see Details.
  }
  \item{Q}{
    a new point pattern dataset or quadrature scheme
    to which the model should be fitted. An object of class
    \code{"ppp"} or \code{"quad"}.
  }
  \item{trend}{
    a new \code{formula} for the spatial trend of the model.
  }
  \item{interaction}{
    a new interpoint interaction structure (object of class
    \code{"interact"}) for the model. 
  }
  \item{covariates}{
    a new data frame of spatial covariates,
    or \code{NULL} to remove all spatial covariates.
  }
  \item{correction}{
    character string giving a new type of edge correction,
    or \code{NULL} to return to the default.
  }
  \item{rbord}{
    numerical value of the distance for the `border' edge correction,
    or \code{NULL} to suppress the border correction.
  }
  \item{use.gam}{
    logical value indicating whether to fit the model using
    \code{\link{gam}}. See \code{\link{ppm}}.
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

  The named arguments \code{Q}, \code{trend},
  \code{interaction}, \code{covariates}, \code{correction},
  \code{rbord}, \code{use.gam} are passed to \code{\link{ppm}}.
  Use \code{name=NULL} to remove the argument \code{name} from the
  call.

  The unnamed arguments \code{...} may include
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
