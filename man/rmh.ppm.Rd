\name{rmh.ppm}
\alias{rmh.ppm}
\title{Simulate from a Fitted Point Process Model}
\description{
  Given a point process model fitted to data, 
  generate a random simulation of the model, 
  using the Metropolis-Hastings algorithm.
}

\usage{rmh.ppm(model,start=NULL,control=NULL,\dots, verbose=TRUE, project=TRUE) }

\arguments{
  \item{model}{A fitted point process model (object of class
    \code{"ppm"}, see \code{\link{ppm.object}}) which it is desired
    to simulate.  This fitted model is usually the result of a call
    to \code{\link{mpl}}.  See \bold{Details} below.
  }
  \item{start}{A list of arguments determining the initial state
    of the Metropolis-Hastings algorithm.  See
    \code{\link{rmh.default}} for description of these arguments.
    Defaults to list(x.start=model$Q$data)
  }
  \item{control}{A list of arguments controlling the running of
    the Metropolis-Hastings algorithm.  See \code{\link{rmh.default}}
    for description of these arguments.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{rmh.default}}.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
  }
  \item{project}{
    Logical flag indicating what to do if the fitted model is
    invalid (in the sense that the values of the fitted coefficients do not
    specify a valid point process).
    If \code{project=TRUE} the closest valid model will be simulated;
    if \code{project=FALSE} an error will occur.
  }
}

\value{A point pattern (an object of class \code{"ppp"}; see
  \code{\link{ppp.object}}).
}

\details{
  This function generates simulated realisations from a point
  process model that has been fitted to point pattern data.  It is
  a method for the generic function \code{\link{rmh}} for the
  class \code{"ppm"} of fitted point process models.  To simulate
  other kinds of point process models, see \code{\link{rmh}}
  or \code{\link{rmh.default}}.

  The argument \code{model} describes the fitted model.  It must be
  an object of class \code{"ppm"} (see \code{\link{ppm.object}}),
  and will typically be the result of a call to the point process
  model fitting function \code{\link{mpl}}.

  The current implementation enables simulation from any fitted model
  involving the interactions
  \code{\link{DiggleGratton}},
  \code{\link{Geyer}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}}
  and \code{\link{Softcore}},
  including nonstationary models, provided the trend does not
  involve spatial covariates. See the examples.

  It is possible that the fitted coefficients of a point process model
  may be ``illegal'', i.e. that there may not exist a
  mathematically well-defined point process with the given parameter
  values. For example, a Strauss process with interaction
  parameter \eqn{\gamma > 1}{gamma > 1} does not exist,
  but the model-fitting procedure used in \code{\link{mpl}} will sometimes
  produce values of \eqn{\gamma}{gamma} greater than 1.
  In such cases, if \code{project=FALSE} then an error will occur,
  while if \code{project=TRUE} then \code{rmh.ppm} will find
  the nearest legal model and simulate
  this model instead. (The nearest legal model is obtained by
  projecting the vector of coefficients onto the set of
  valid coefficient vectors. The result is usually the Poisson process
  with the same fitted intensity.)
  
  The arguments \code{start} and \code{control} are lists of
  parameters determining the initial state and the iterative
  behaviour, respectively, of the Metropolis-Hastings algorithm.
  They are passed directly to \code{\link{rmh.default}}.
  See \code{\link{rmh.default}} for details of these parameters.
  Note that if the model has a trend then the \code{expand}
  component of \code{control} will be ignored.
  The value of \code{expand} must be equal to 1 in this setting.

  After extracting the relevant information from the fitted model
  object \code{model}, \code{rmh.ppm} simply invokes the default
  \code{rmh} algorithm \code{\link{rmh.default}}, unless the model
  is Poisson.  If the model is Poisson then \code{rmh.default} is
  not needed, and a stand-in function \code{pseudo.rmh} is invoked
  to do the simple simulation that is required.

  See \code{\link{rmh.default}} for further information about the
  implementation, or about the Metropolis-Hastings algorithm.
}

\section{Warnings}{
  See Warnings in \code{\link{rmh.default}}.
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
  \code{\link{Geyer}},
  \code{\link{DiggleGratton}}
}

\examples{
   require(spatstat)
   
   data(swedishpines)
   X <- swedishpines
   plot(X, main="Swedish Pines data")

   # Poisson process
   fit <- mpl(X, ~1, Poisson())
   Xsim <- rmh(fit)
   plot(Xsim, main="simulation from fitted Poisson model")

   # Strauss process   
   fit <- mpl(X, ~1, Strauss(r=7), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Strauss model")

   # Strauss - hard core process
   fit <- mpl(X, ~1, StraussHard(r=7,hc=2), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Strauss hard core model")

   # Geyer saturation process
   fit <- mpl(X, ~1, Geyer(r=7,sat=2), rbord=7)
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Geyer model")

   # soft core interaction process
   fit <- mpl(X, ~1, Softcore(kappa=0.1))
   Xsim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Soft Core model")

   data(cells)
   plot(cells)
   # Diggle-Gratton pairwise interaction model
   fit <- mpl(cells, ~1, DiggleGratton(0.05, 0.1))
   Xsim <- rmh(fit, start=list(n.start=cells$n), control=list(nrep=1e3))
   plot(Xsim, main="simulation from fitted Diggle-Gratton model")
   
   # marked point pattern
   data(amacrine)
   Y <- amacrine
   plot(Y, main="Amacrine data")

   # marked Poisson models 
   fit <- mpl(Y)
   Ysim <- rmh(fit)
   plot(Ysim, main="simulation from mpl(Y)")

   fit <- mpl(Y,~marks)
   Ysim <- rmh(fit)
   plot(Ysim, main="simulation from mpl(Y, ~marks)")

   fit <- mpl(Y,~polynom(x,y,2))
   Ysim <- rmh(fit)
   plot(Ysim, main="simulation from mpl(Y, ~polynom(x,y,2))")

   fit <- mpl(Y,~marks+polynom(x,y,2))
   Ysim <- rmh(fit)
   plot(Ysim, main="simulation from mpl(Y, ~marks+polynom(x,y,2))")

   fit <- mpl(Y,~marks*polynom(x,y,2))
   Ysim <- rmh(fit)
   plot(Ysim, main="simulation from mpl(Y, ~marks*polynom(x,y,2))")

   # multitype Strauss models
   MS <- MultiStrauss(types = levels(Y$marks),
                      radii=matrix(0.07, ncol=2, nrow=2))
   fit <- mpl(Y, ~marks, MS)
   Ysim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Ysim, main="simulation from fitted Multitype Strauss")

   fit <- mpl(Y,~marks*polynom(x,y,2), MS)
   Ysim <- rmh(fit, start=list(n.start=X$n), control=list(nrep=1e3))
   plot(Ysim, main="simulation from fitted inhomogeneous Multitype Strauss")

}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
