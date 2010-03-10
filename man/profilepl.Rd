\name{profilepl}
\alias{profilepl}
\title{Profile Maximum Pseudolikelihood}
\description{
  Fits point process models by profile maximum pseudolikelihood
}
\usage{
profilepl(s, f, \dots, rbord = NULL, verbose = TRUE)
}
\arguments{
  \item{s}{
    Data frame containing values of the irregular parameters
    over which the profile pseudolikelihood will be computed.
  }
  \item{f}{
    Function (such as \code{\link{Strauss}})
    that generates an interpoint interaction object, given
    the values of the irregular parameters.
  }
  \item{\dots}{
    Data passed to \code{\link{ppm}} to fit the model.
  }
  \item{rbord}{
    Radius for border correction (same for all models).
    If omitted, this will be computed from the interactions.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
}
\details{
  The model-fitting function \code{\link{ppm}} fits point process
  models to point pattern data. However, 
  only the \sQuote{regular} parameters of the model can be fitted by
  \code{\link{ppm}}. The model may also depend on \sQuote{irregular}
  parameters that must be fixed in any call to \code{\link{ppm}}.

  This function \code{profilepl} is a wrapper which finds the values of the
  irregular parameters that give the best fit. It uses the method of
  maximum profile pseudolikelihood. 

  The argument \code{f} would typically be one of the functions
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{DiggleGratton}},
  \code{\link{Geyer}},
  \code{\link{LennardJones}}
  or \code{\link{OrdThresh}}.
  For the moment, assume this is so.
  
  The argument \code{s} must be a data frame whose columns contain
  values of the irregular parameters. The names of the columns of
  \code{s} must match the argument names of \code{f}.

  To apply the method of profile maximum pseudolikelihood,
  each row of \code{s} will be taken in turn; the parameter values in this row
  will be passed to \code{f}, resulting in an interaction object.
  Then \code{\link{ppm}} will be applied to the data \code{...}
  using this interaction; this results in a fitted point process model.
  The value of the log pseudolikelihood from this model is stored.
  After all rows of \code{s} have been processed in this way, the
  row giving the maximum value of log pseudolikelihood will be found.

  The object returned by \code{profilepl} contains the profile
  pseudolikelihood function, the best fitting model, and other data.
  It can be plotted (yielding a
  plot of the log pseudolikelihood values against the irregular
  parameters) or printed (yielding information about the best fitting
  values of the irregular parameters). 
 
  In general, \code{f} may be any function that will return
  an interaction object (object of class \code{"interact"})
  that can be used in a call to \code{\link{ppm}}. Each argument of
  \code{f} must be a single value. 
}
\value{
  An object of class \code{"profilepl"}. There are methods
  for \code{\link{plot}} and \code{\link{print}} for this class.

  The components of the object include
  \item{fit}{Best-fitting model}
  \item{param}{The data frame \code{s}}
}
\examples{
    data(cells)

    # one irregular parameter
    s <- data.frame(r=seq(0.05,0.15, by=0.01))
    \testonly{
      s <- data.frame(r=c(0.05,0.1,0.15))
    }
    ps <- profilepl(s, Strauss, cells)
    ps
    plot(ps)

    # two irregular parameters
    s <- expand.grid(r=seq(0.05,0.15, by=0.01),sat=1:3)
    \testonly{
      s <- expand.grid(r=c(0.05,0.1,0.15),sat=1:2)
    }
    pg <- profilepl(s, Geyer, cells)
    pg
    plot(pg)
    \dontrun{
    pg$fit
    }
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}

