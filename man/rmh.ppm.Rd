\name{rmh.ppm}
\alias{rmh.ppm}
\title{Simulate from a Fitted Point Process Model}
\description{
  Given a point process model fitted to data, 
  generate a random simulation of the model, 
  using the Metropolis-Hastings algorithm.
}

\usage{rmh.ppm(model,start,control,\dots) }

\arguments{
  \item{model}{A fitted point process model
    (object of class \code{"ppm"}, see \code{\link{ppm.object}})
    which it is desired to simulate.
    This fitted model is usually the result of a call to \code{\link{mpl}}.
    See \bold{Details} below.
  }
  \item{start}{A list of arguments determining the initial
    state of the Metropolis-Hastings algorithm.
    See \code{\link{rmh.default}} for description of these arguments.
  }
  \item{control}{A list of arguments controlling the
    running of the Metropolis-Hastings algorithm.
    See \code{\link{rmh.default}} for description of these arguments.
  }
  \item{\dots}{Further arguments are ignored.}
}

\value{A point pattern (an object of class \code{"ppp"}, see
  \code{\link{ppp.object}}).
}

\details{
  This function generates simulated realisations
  from a point process model that has been fitted to point pattern data.
  It is a method for the generic function \code{\link{rmh}}
  for the class \code{"ppm"} of fitted point process models.
  To simulate other kinds of point process models,
  see \code{\link{rmh}} or \code{\link{rmh.default}}.
  
  The argument \code{model} describes the fitted model.
  It must be an object of class \code{"ppm"} (see
  \code{\link{ppm.object}}) and will typically 
  be the result of a call to the point process model fitting
  function \code{\link{mpl}}.

  The current implementation is experimental, and
  only a few processes can be simulated.
  At present the fitted model must not have any spatial trend,
  and the only models possible are 
  the Poisson, Strauss, Strauss/Hard Core, Soft Core, and Geyer interactions.
  These are fitted by \code{\link{mpl}} using
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}}
  and \code{\link{Geyer}} respectively.
  See the examples.

  The arguments \code{start} and \code{control}
  are lists of parameters determining
  the initial state and the iterative behaviour, respectively,
  of the Metropolis-Hastings algorithm.
  They are passed directly to \code{\link{rmh.default}}.
  See \code{\link{rmh.default}} for details of these parameters.

  After extracting the relevant
  information from the fitted model object \code{model},
  \code{rmh.ppm} simply invokes the default \code{rmh} algorithm
  \code{\link{rmh.default}}.

  See \code{\link{rmh.default}} for further information
  about the implementation, or about the Metropolis-Hastings algorithm.
}

\section{Warnings}{
See Warnings in \code{\link{rmh.default}}
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmh.default}},
  \code{\link{ppp.object}},
  \code{\link{mpl}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{Geyer}}
}

\examples{
   require(spatstat)
   data(swedishpines)
   X <- swedishpines
   plot(X, main="Swedish Pines data")

   fit <- mpl(X, ~1, Strauss(r=7), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Strauss model")

   fit <- mpl(X, ~1, StraussHard(r=7,hc=2), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Strauss hard core model")

   fit <- mpl(X, ~1, Geyer(r=7,sat=2), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Geyer model")

   fit <- mpl(X, ~1, Softcore(kappa=0.1))
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Soft Core model")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
