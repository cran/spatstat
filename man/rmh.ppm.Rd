\name{rmh.ppm}
\alias{rmh.ppm}
\title{Simulate from a Fitted Point Process Model}
\description{
  Given a point process model fitted to data, 
  generate a random simulation of the model, 
  using the Metropolis-Hastings algorithm.
}

\usage{
  \method{rmh}{ppm}(model, start=NULL,
                    control=default.rmhcontrol(model, expand=expand),
                    \dots,
                    project=TRUE, expand=NULL, 
                    verbose=TRUE)
}

\arguments{
  \item{model}{A fitted point process model (object of class
    \code{"ppm"}, see \code{\link{ppm.object}}) which it is desired
    to simulate.  This fitted model is usually the result of a call
    to \code{\link{ppm}}.  See \bold{Details} below.
  }
  \item{start}{Data determining the initial state
    of the Metropolis-Hastings algorithm.  See
    \code{\link{rmhstart}} for description of these arguments.
    Defaults to \code{list(x.start=data.ppm(model))}
  }
  \item{control}{Data controlling the running of
    the Metropolis-Hastings algorithm.  See \code{\link{rmhcontrol}}
    for description of these arguments.
  }
  \item{\dots}{
    Further arguments passed to \code{\link{rmh.default}}.
  }
  \item{project}{
    Logical flag indicating what to do if the fitted model is
    invalid (in the sense that the values of the fitted coefficients do not
    specify a valid point process).
    If \code{project=TRUE} the closest valid model will be simulated;
    if \code{project=FALSE} an error will occur.
  }
  \item{expand}{
    Optional expansion factor for simulation window.
    Defaults to 2.
  }
  \item{verbose}{
    Logical flag indicating whether to print progress reports.
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
  model fitting function \code{\link{ppm}}.

  The current implementation enables simulation from any fitted model
  involving the interactions
  \code{\link{DiggleGratton}},
  \code{\link{Geyer}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}}
  and \code{\link{Softcore}},
  including nonstationary models. See the examples.

  It is possible that the fitted coefficients of a point process model
  may be ``illegal'', i.e. that there may not exist a
  mathematically well-defined point process with the given parameter
  values. For example, a Strauss process with interaction
  parameter \eqn{\gamma > 1}{gamma > 1} does not exist,
  but the model-fitting procedure used in \code{\link{ppm}} will sometimes
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
  They are passed directly to \code{\link{rmhstart}} and
  \code{\link{rmhcontrol}} respectively.
  See \code{\link{rmhstart}} and
  \code{\link{rmhcontrol}} for details of these parameters.

  Note that if you specify \code{control$expand > 1} (so that the
  model will be simulated on a window larger than the original data
  window) then the model must be capable of extrapolation to this
  larger window. This excludes models which depend on external covariates.
  
  After extracting the relevant information from the fitted model
  object \code{model}, \code{rmh.ppm} simply invokes the default
  \code{rmh} algorithm \code{\link{rmh.default}}, unless the model
  is Poisson.

  If the model is Poisson then the Metropolis-Hastings
  algorithm is not needed, and the model is simulated directly, using
  one of \code{\link{rpoispp}}, \code{\link{rmpoispp}},
  \code{\link{rpoint}} or \code{\link{rmpoint}}.

  See \code{\link{rmh.default}} for further information about the
  implementation, or about the Metropolis-Hastings algorithm.
}

\section{Warnings}{
  See Warnings in \code{\link{rmh.default}}.
}

\seealso{
  \code{\link{simulate.ppm}},
  \code{\link{rmh}},
  \code{\link{rmhmodel}},
  \code{\link{rmhcontrol}},
  \code{\link{default.rmhcontrol}},
  \code{\link{rmhstart}},
  \code{\link{rmh.default}},
  \code{\link{ppp.object}},
  \code{\link{ppm}},
  \code{\link{PairPiece}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{Geyer}},
  \code{\link{AreaInter}},
  \code{\link{DiggleGratton}}
}

\examples{

   op <- spatstat.options(rmh.nrep=1e5)
   live <- interactive()

   \testonly{
     spatstat.options(rmh.nrep=10)
   }

   data(swedishpines)
   X <- swedishpines
   plot(X, main="Swedish Pines data")

   # Poisson process
   fit <- ppm(X, ~1, Poisson())
   Xsim <- rmh(fit)
   if(live) plot(Xsim, main="simulation from fitted Poisson model")

   # Strauss process   
   fit <- ppm(X, ~1, Strauss(r=7))
   Xsim <- rmh(fit)
   if(live) plot(Xsim, main="simulation from fitted Strauss model")

   # Strauss process simulated on a larger window
   # then clipped to original window
   Xsim <- rmh(fit, control=list(nrep=1e3, expand=2, periodic=TRUE))
   
   # Strauss - hard core process
   fit <- ppm(X, ~1, StraussHard(r=7,hc=2))
   Xsim <- rmh(fit, start=list(n.start=X$n))
   if(live) plot(Xsim, main="simulation from fitted Strauss hard core model")

   # Geyer saturation process
   fit <- ppm(X, ~1, Geyer(r=7,sat=2))
   Xsim <- rmh(fit, start=list(n.start=X$n))
   if(live) plot(Xsim, main="simulation from fitted Geyer model")

   # Area-interaction process
   fit <- ppm(X, ~1, AreaInter(r=7))
   Xsim <- rmh(fit, start=list(n.start=X$n))
   if(live) plot(Xsim, main="simulation from fitted area-interaction model")

   # soft core interaction process
   Q <- quadscheme(X, nd=50)
   fit <- ppm(Q, ~1, Softcore(kappa=0.1), correction="isotropic")
   Xsim <- rmh(fit, start=list(n.start=X$n))
   if(live) plot(Xsim, main="simulation from fitted Soft Core model")

   data(cells)
   if(live) plot(cells)
   # Diggle-Gratton pairwise interaction model
   fit <- ppm(cells, ~1, DiggleGratton(0.05, 0.1))
   Xsim <- rmh(fit, start=list(n.start=cells$n))
   if(live) plot(Xsim, main="simulation from fitted Diggle-Gratton model")

   X <- rSSI(0.05, 100)
   if(live) plot(X, main="new data")

   # piecewise-constant pairwise interaction function
   fit <- ppm(X, ~1, PairPiece(seq(0.02, 0.1, by=0.01)))
   Xsim <- rmh(fit)
   if(live) plot(Xsim, main="simulation from fitted pairwise model")

   # marked point pattern
   data(amacrine)
   Y <- amacrine
   if(live) plot(Y, main="Amacrine data")

   # marked Poisson models 
   fit <- ppm(Y)
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from ppm(Y)")

   fit <- ppm(Y,~marks)
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from ppm(Y, ~marks)")

   fit <- ppm(Y,~polynom(x,y,2))
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from ppm(Y, ~polynom(x,y,2))")

   fit <- ppm(Y,~marks+polynom(x,y,2))
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from ppm(Y, ~marks+polynom(x,y,2))")

   fit <- ppm(Y,~marks*polynom(x,y,2))
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from ppm(Y, ~marks*polynom(x,y,2))")

   # multitype Strauss models
   MS <- MultiStrauss(types = levels(Y$marks),
                      radii=matrix(0.07, ncol=2, nrow=2))
   fit <- ppm(Y, ~marks, MS)
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from fitted Multitype Strauss")

   fit <- ppm(Y,~marks*polynom(x,y,2), MS)
   Ysim <- rmh(fit)
   if(live) plot(Ysim, main="simulation from fitted inhomogeneous Multitype Strauss")

   spatstat.options(op)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{models}
\keyword{datagen}
