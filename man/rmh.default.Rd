\name{rmh.default}
\alias{rmh.default}
\title{Simulate Point Process Models using the Metropolis-Hastings Algorithm.}
\description{
  Generates a random point pattern, simulated from
  a chosen point process model, using the Metropolis-Hastings
  algorithm. 
}

\usage{rmh.default(model,start,control,verbose=TRUE,\dots) }

\arguments{
  \item{model}{A named list of objects specifying the point process model
    that is to be simulated, having (some of) the following components:
    \describe{
      \item{cif}{Character string specifying the choice of model}
      \item{par}{Parameters of the model}
      \item{w}{Spatial window in which to simulate}
      \item{trend}{Either a \bold{function} whose first two arguments
      are \code{x} and \code{y} (the coordinates of an abitrary
      point in the observation window) which calculates the value
      of the trend at the point \code{(x,y)} or an \bold{image}
      which specifies the value of the trend at (a grid of) points
      in the observation window.  Auxilliary arguments may be passed
      to the \code{trend} function via the \dots argument.
      If the model being simulated is multitype, then \code{trend}
      may be a \bold{list} of functions or images, with as many entries
      as there are types.}
      \item{types}{A vector of the possible levels of a factor forming
      the marks for a multitype process.  If the model being simulated
      is multitype and \code{types} is not specified, then this vector
      defaults to \code{1:ntypes} where \code{ntypes} is the number of
      types.}
    }
    See \emph{Details} for details.
  }
  \item{start}{List of parameters determining the initial state of
    the algorithm:
    \describe{
      \item{n.start}{number of initial points (randomly generated)}
      \item{x.start}{initial point pattern configuration}
      \item{iseed}{vector of seeds for random number generator}
    }
    The parameters \code{n.start} and \code{x.start} are incompatible;
    precisely \bold{one} of them should be specified.
    See \emph{Details} for details.
  }
  \item{control}{List of parameters controlling the iterative behaviour
    and termination of the algorithm:
    \describe{
      \item{p}{Probability of proposing a shift (as against a birth/death)}
      \item{q}{Conditional probability of proposing a death given that a
	birth or death will be proposed}
      \item{nrep}{Total number of steps (proposals) of Metropolis-Hastings
	algorithm that should be run}
      \item{expand}{Factor by which the window will be expanded}
      \item{periodic}{(Logical) whether to apply periodic boundary
	conditions}
      \item{ptypes}{For multitype point processes, the distribution of
	the mark attached to a new random point (when a birth is
	proposed)}
      \item{fixall}{(Logical) for multitype point processes,
	whether to fix the number of points of each type.}
      \item{nverb}{Progress reports will be printed every \code{nverb}
	iterations}
    }
    See \emph{Details} for details.
  }
  \item{verbose}{A logical scalar; if \code{verbose} is TRUE then
   warnings are printed out whenever the storage space alloted to
   the underlying Fortran code, to hold the generated points, gets
   increased.  If it is FALSE then this process proceeds silently.}
  \item{\dots}{Further arguments, e.g. to be passed to \bold{trend}
    functions.}
}

\value{A point pattern (an object of class \code{"ppp"}, see
  \code{\link{ppp.object}}).

  The returned value has an attribute \code{info} consisting of
  arguments supplied to the function (or default values of arguments
  which were not explicitly supplied).  These are given so that it
  is possible to reconstruct exactly the manner in which the pattern
  was generated.  The components of \code{info} are \code{model},
  \code{start}, and \code{control} which in turn are lists:

  \code{model=list(cif, par, trend)}
  
  \code{start=list(n.start,x.start,iseed)}

  \code{control=list(p=p,q=q,nrep=nrep,expand,periodic,
                     ptypes=ptypes,fixall=fixall)}

  Note that only one of \code{x.start} and \code{x.start} appear in
  in the \code{start} list.
}

\details{
  This function generates simulated realisations from any of a range of
  spatial point processes, using the Metropolis-Hastings algorithm.
  It is the default method for the generic function \code{rmh}.

  This function executes a Metropolis-Hastings algorithm
  with birth, death and shift proposals as described in
  Geyer and Moller (1994).

  The argument \code{model} specifies the point process model to be
  simulated. It is a list with the following components:

  \describe{
    \item{cif}{A character string specifying the choice of
      interpoint interaction for the point process. The current options are
      \describe{
	\item{\code{'strauss'}}{The Strauss process}
	\item{\code{'straush'}}{The Strauss process with hard core}
	\item{\code{'sftcr'}}{The Softcore process}
	\item{\code{'straussm'}}{ The multitype Strauss process}
        \item{\code{'straushm'}}{Multitype Strauss process with hard core}
        \item{\code{'dgs'}}{Diggle, Gates and Stibbard (1987) process}
	\item{\code{'diggra'}}{Diggle and Gratton (1984) process}
	\item{\code{'geyer'}}{Saturation process (Geyer, 1999).}
      }
      See the section \bold{Extensions} for the possibility of extending
      this list of options.
    }
    \item{par}{
      A vector (or a list if the pattern is multitype) providing
      a set of parameter values appropriate to the conditional
      intensity function being invoked.  The structure of \code{par}
      should be as follows, for each of the available conditional
      intensity functions:

      \describe{
	\item{strauss:}{
	  (Strauss process.) A \bold{named} vector with components
	  \code{beta,gamma,r} which are respectively the ``base''
	  intensity, the pair-wise interaction parameter and the
	  interaction radius.  Note that \code{gamma} must be less than
	  or equal to 1.
	}
	\item{straush:}{
	  (Strauss process with hardcore.) A \bold{named} vector with
	  entries \code{beta,gamma,r,hc} where \code{beta}, \code{gamma},
	  and \code{r} are as for the Strauss process, and \code{hc} is
	  the hardcore radius.  Of course \code{hc} must be less than
	  \code{r}.
	}
	\item{sftcr:}{
	  (Softcore process.) A \bold{named} vector with components
	  \code{beta,sigma,kappa}.  Again
	  \code{beta} is a ``base'' intensity.  The pairwise
	  interaction between two points \eqn{u \neq v}{u != v} is
	  \deqn{
	    \exp \left \{
	    - \left ( \frac{\sigma}{||u-v||} \right )^{2/\kappa}
	    \right \}
	  }{-(sigma/||u-v||)^(2/kappa)}
	  Note that it is necessary that \eqn{0 < \kappa < 1}{0 < kappa <1}.
	}
	\item{straussm:}{
	  (Multitype Strauss process.) A \bold{named} list with components
	  \itemize{
	    \item
	    \code{beta}: 
	      A vector of ``base'' intensities, one for each possible
	      type.
	    \item
	    \code{gamma}:
	      A \bold{symmetric} matrix of interaction parameters,
	      with \eqn{\gamma_{ij}}{gamma_ij} pertaining to the interaction between
	      type \eqn{i} and type \eqn{j}.
	    \item
	    \code{radii}:
	      A \bold{symmetric} matrix of interaction radii, with
	      entries \eqn{r_{ij}}{r_ij} pertaining to the interaction between type
	      \eqn{i} and type \eqn{j}.
	   }
	 }
	\item{straushm:}{
	  (Multitype Strauss process with hardcore.)
	  A \bold{named} list with components \code{beta} and \code{gamma}
	  as for \code{straussm} and
	  \bold{two} ``radii'' components:
	  \itemize{
            \item \code{iradii}: the interaction radii
            \item \code{hradii}: the hardcore radii
	  }
	  which are both symmetric matrices of nonnegative numbers.
	  The entries of \code{hradii} must be less than the
	  corresponding entries
	  of \code{iradii}.
	}
	\item{dgs:}{
	  (Diggle, Gates, and Stibbard process.
	  See Diggle, Gates, and Stibbard (1987))
	  A \bold{named} vector with components
	  \code{beta} and \code{rho}.  This process has pairwise interaction
	  function equal to
	  \deqn{
	    e(t) = \sin^2\left(\frac{\pi t}{2\rho}\right)
	  }{
	    e(t) = sin^2((pi * t)/(2 * rho))
	  }
	  for \eqn{t < \rho}{t < rho}, and equal to 1
	  for \eqn{t \ge \rho}{t >= rho}.
	}
	\item{diggra:}{
	  (Diggle-Gratton process. See Diggle and Gratton (1984)
	  and Diggle, Gates and Stibbard (1987).)
	  A \bold{named} vector with components \code{beta},
	  \code{kappa}, \code{delta} and \code{rho}.  This process has
	  pairwise interaction function \eqn{e(t)} equal to 0
	  for \eqn{t < \delta}{t < delta}, equal to
	  \deqn{
	    \left(\frac{t-\delta}{\rho-\delta}\right)^\kappa
	  }{
	    ((t-delta)/(rho-delta))^kappa
	  }
	  for \eqn{\delta \le t < \rho}{delta <= t < rho},
	  and equal to 1 for \eqn{\delta \ge  \rho}{delta >= rho}.
	  Note that here we use the symbol
	  \eqn{\kappa}{kappa} where Diggle, Gates, and Stibbard use
	  \eqn{\beta}{beta} since we reserve the symbol \eqn{\beta}{beta}
	  for an intensity parameter.
	}
	\item{geyer}{
	  (See Geyer (1999))  A \bold{named} vector
	  with components \code{beta}, \code{gamma}, \code{r}, and \code{sat}.
	  The components \code{beta}, \code{gamma}, \code{r} are as for
	  the Strauss model, and \code{sat} is the ``saturation''
	  parameter.  The model is Geyer's ``saturation'' point process
	  model.  It is ``like a Strauss model, but with an upper bound
	  to the number of \eqn{r}-close neighbors of any point.''
	 
	  Explicitly, a saturation point process with interaction
	  radius \eqn{r}, saturation threshold \eqn{s}, and  parameters
	  \eqn{\beta}{beta} and \eqn{\gamma}{gamma}, is the point process
	  in which each point \eqn{x_i}{x[i]}
	  in the pattern \eqn{X} contributes a factor
	  \deqn{\beta \gamma^{\min(s, t(x_i,X))}}{beta gamma^min(s,t(x[i],X))}
	  to the probability density of the point pattern, where
	  \eqn{t(x_i,X)}{t(x[i],X)}
	  denotes the number of ``\eqn{r}-close neighbours''
	  of \eqn{x_i}{x[i]} in the pattern \eqn{X}.
	}
      }
    }
    \item{trend}{
	A function or image (or a list of such if the model
	is multitype) to provide the value of the trend
	(if there is one) at an arbitrary point.  A trend
	function may be a function of any number of arguments,
	but the first two must be the spatial coordinates of
	a point.  The function \bold{must} be \bold{vectorized}.
	That is, it must be capable of accepting vector valued
	\code{x} and \code{y} arguments.  Put another way,
	it must be capable of calculating the trend value at a
	number of points, simultaneously, and should return the
	\bold{vector} of corresponding trend values.  An image (see
	\code{\link{im}()}) provides the trend values at a grid of
	points in the observation window and determines the trend
	value at other points as the value at the nearest grid point.

	Note that the trend or trends must be \bold{non-negative};
	no checking is done for this.

    }
    \item{w}{
      A specification of a window in which the pattern is
      to be generated.  If specified, it must be in a form
      which can be coerced to an object of class \code{owin}
      by \code{\link{as.owin}()}.
    }
  }
  The argument \code{start} determines the initial state of the
  Metropolis-Hastings algorithm. Possible components are
    \describe{
      \item{n.start}{The number of ``initial'' points to be randomly
	(uniformly) generated in the window \code{w}.  For a
	multitype point process, \code{n.start} may be a vector
	(of length equal to the number of types) giving the number
	of points of each type to be generated.  Incompatible with
	\code{x.start}.

	The component \code{n.start} may be a scalar (integer)
        or, for a multitype process, a vector of integers of
	length equal to the number of types.

	A vector-valued \code{n.start} is meaningful \bold{only}
	if \code{p} (the probability of a shift as opposed to a birth or
	death) is equal to 1 (so that we are conditioning on the number
	of points).  When \eqn{p < 1}, if \code{n.start} is vector valued
	then it is effectively replaced by its sum.

	The resulting set of uniformly generated points gives
	the Metropolis-Hastings algorithm an initial state from
	which to start.  (Actually, when \eqn{p < 1}, the number
	\code{n.start} gets multiplied by the ratio of the area of
	the \bold{expanded} window to that of the original window.
	Then that many points are uniformly generated in the expanded
	window; see below for a discussion of the expanded window.
	The value of \code{n.start} should be roughly equal to
	(an educated guess at) the expected number of points which
	will be generated inside the window.
      }
      \item{x.start}{
	Initial point pattern configuration.  A point pattern (an
	object of class \code{"ppp"}, or data which can be coerced
	to this class by \code{\link{as.ppp}}).  Incompatible with
	\code{n.start}.

	The component \code{x.start} is a point pattern (an object
	of class \code{ppp}, or an object which can be coerced
	to this class by \code{as.ppp()}).  This object provides
	an alternative means of specifying the initial ``state''
	or configuration for the Metropolis-Hastings algorithm.
	If \code{x.start} is specified, but \code{x.start$window}
	is \code{NULL}, then this gap is filled in by the component
	\code{w} of \code{start}.  If \code{x.start$window}
	is present, and if \code{w} is specified as well, the
	latter is used as a window to which to clip the final
	simulated pattern.  Thus in such circumstances it is only
	sensible to specify a value of \code{w} which is contained
	in \code{x.start$window}.  However \bold{no checking}
	is done for this.

	The simulated pattern is constructed in \code{x.start$window}.
	No expansion takes place.  (The argument \code{expand} is
	forced to equal 1.)  As indicated above, at the end of the
	simulation, the resulting pattern is clipped to the window
	\code{w} if this is given.
      }
      \item{iseed}{Seed for random number generator.
	A triple of integers.  If unspecified these are themselves
	generated, on the interval from 1 to 1 million, using the
	function \code{\link{sample}()}.
      }
    }
    The parameters \code{n.start} and \code{x.start} are \emph{incompatible}.

    The third argument \code{control} controls the simulation
    procedure, iterative behaviour, and termination of the
    Metropolis-Hastings algorithm. It is a list with components:
     \describe{
       \item{p}{The probability of of proposing a ``shift''
	 (as opposed to a birth or death) in the Metropolis-Hastings
	 algorithm.  If \eqn{p = 1} then we do nothing but shifts,
	 whence the number of points never changes, i.e. we are
	 simulating conditionally upon the number of points.  In this
	 case, for multitype processes, we also have the option of
	 simulating conditionally upon the number of points of each
	 type; this is effected by specifying \code{fixall} (see
	 below) to equal \code{TRUE}.  In this case, \code{n.start}
	 must be a vector whose entries are these numbers.

	 We can only condition on the number of points if the simulation
	 takes place in the original window (as opposed to taking place in
	 a larger window and then being clipped to the original).  Hence,
	 if \eqn{p = 1} then \code{expand} defaults to 1 and it is an
	 error to specify a value of \code{expand} which is greater than 1.
       }
       \item{q}{The probability of proposing a death (rather than a birth)
	 given that birth/death has been chosen over shift.  This is of
	 course ignored if \code{p} is equal to 1.
       }
       \item{nrep}{The number of repetitions or steps (changes of state)
	 to be made by the Metropolis-Hastings algorithm.  It should
	 be large.
       }
       \item{expand}{The factor by which the enclosing box of the window
	 \code{w} is to be expanded in order to better approximate
	 the simulation of a process existing in the whole plane,
	 rather than just in the window.  If \code{expand} equals 1,
	 then we are simulating the latter, unless \code{periodic}
	 (see below) is \code{TRUE}.  The larger \code{expand} is,
	 the better we approximate the former.	Note that any value
	 of \code{expand} smaller than 1 is treated as if it were 1.

	 The area of the expanded window is equal to \code{expand}
	 times the area of the enclosing box; width and height are
	 stretched proportionately.  Points are generated by the
	 Metropolis-Hastings algorithm in the expanded window, and
	 then ``clipped'' down to the original window when algorithm
	 has finished.	The argument \code{expand} defaults to 2
	 if \code{periodic} is \code{FALSE} and \eqn{p < 1} and
	 to 1 if \code{periodic} is \code{TRUE} or if \eqn{p <
	 1}, or if the starting configuration is specified via
	 \code{x.start}.  Trying to set \code{expand} greater
	 than 1 when \code{periodic} is \code{TRUE} or \eqn{p = 1}
	 generates an error.  A specified value of \code{expand}
	 is simply ignored if \code{x.start} is given.
       }
       \item{periodic}{A logical scalar; if \code{periodic} is \code{TRUE}
	 we simulate a process on the torus formed by
	 identifying opposite edges of the (rectangular) window.
	 If \code{periodic} is \code{TRUE} and the window is not
	 rectangular, an error is given.
       }
       \item{ptypes}{A vector of probabilities (summing to 1) to be used
	 in assigning a random type to a new point.  Defaults to a vector
	 each of whose entries is \eqn{1/nt} where \eqn{nt} is the number
	 of types for the process.  Convergence of the simulation
	 algorithm should be improved if \code{ptypes} is close to the
	 relative frequencies of the types which will result from the
	 simulation.
       }
       \item{fixall}{A logical scalar specifying whether to condition on
	 the number of points of each type.  Meaningful only if a marked
	 process is being simulated, and if \eqn{p = 1}.  A warning message
	 is given if \code{fixall} is set equal to \code{TRUE} when it is
	 not meaningful.
       }
       \item{nverb}{An integer specifying how often ``progress reports''
	 (which consist simply of the number of repetitions completed)
	 should be printed out.  If nverb is left at 0, the default,
	 the simulation proceeds silently.
       }
     }
   }
 
\references{
   Baddeley, A. and Turner, R. (2000) Practical maximum
   pseudolikelihood for spatial point patterns.
   \emph{Australian and New Zealand Journal of Statistics}
   \bold{42}, 283 -- 322.

   Diggle, P.J. and Gratton, R.J. (1984)
   Monte Carlo methods of inference for implicit statistical models.
   \emph{Journal of the Royal Statistical Society, series B}
   \bold{46}, 193 -- 212.

   Diggle, P.J., Gates, D.J., and Stibbard, A. (1987)
   A nonparametric estimator for pairwise-interaction point processes.
   Biometrika \bold{74}, 763 -- 770.

   Geyer, C.J. and M{\o}ller, J. (1994)
   Simulation procedures and likelihood inference for spatial
   point processes.
   \emph{Scandinavian Journal of Statistics} \bold{21}, 359--373.

   Geyer, C.J. (1999)
   Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) \emph{Stochastic Geometry: Likelihood and
   Computation}, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79--140.
}

\section{Note}{
  It is possible to simulate conditionally upon the number of
  points, or in the case of multitype processes, upon the number of
  points of each type.  To condition upon the total number of points,
  set \code{p} (the probability of a shift) equal to 1, and specify
  \code{n.start} to be a scalar (as usual).  To condition upon the
  number of points of each type, set \code{p} equal to 1, \code{fixall}
  equal to \code{TRUE}, and specify \code{n.start} to be a vector of
  length \eqn{nt} where \eqn{nt} is the number of types.

  In these circumstances
  \itemize{
    \item
    The value of \code{expand} must be equal to 1; it
    defaults to 1, and it is an error to specify a value larger
    than 1.
    \item
    The resulting simulated pattern will have precisely
    the number of points (of each type) specified by
    \code{n.start}.
  }
}
    
\section{Warnings}{
There is never a guarantee that the Metropolis-Hastings algorithm
has converged to the steady state.

If \code{x.start} is specified then \code{expand} is set equal to 1
and simulation takes place in \code{x.start$window}.  Any specified
value for \code{expand} is simply ignored.

The presence of both a component \code{w} of \code{model} and a
non-null value for \code{x.start$window} makes sense ONLY if \code{w}
is contained in \code{x.start$window}.  However \bold{no checking}
is done for this.

For multitype processes make sure that, even if there is to be no
trend corresponding to a particular type, there is still a component
(a NULL component) for that type, in the list.

No checking is done to verify non-negativity of any specified trend
or trends.
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmh.ppm}},
  \code{\link{ppp}},
  \code{\link{mpl}},
  \code{\link{Strauss}},
  \code{\link{Softcore}},
  \code{\link{StraussHard}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}},
  \code{\link{DiggleGratton}}
}
\section{Extensions}{
  The argument \code{model$cif} matches the name of a Fortran
  subroutine which calculates the conditional intensity function
  for the model.  It is intended that more options will be
  added in the future.  The very brave user \bold{could} try
  to add her own.  Note that in addition to writing Fortran
  code for the new conditional intensity function, the user
  would have to modify the code in the files \code{cif.f} and
  \code{rmh.default.R} appropriately.  (And of course re-install
  the \code{spatstat} package so as to update the dynamically
  loadable shared object \code{spatstat.so}.)
}
\examples{
   \dontrun{
   require(spatstat)
   nr   <- 1e5
   nv  <- 5000
   }
   \testonly{
   nr  <- 10
   nv <- 5
   }
   set.seed(961018)
   
   # Strauss process.
   mod01 <- list(cif="strauss",par=c(beta=2,gamma=0.2,r=0.7),
                 w=c(0,10,0,10))
   X1.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process, conditioning on n = 80:
   X2.strauss <- rmh(model=mod01,start=list(n.start=80),
                     control=list(p=1,nrep=nr,nverb=nv))
   
   # Strauss process equal to pure hardcore:
   mod02 <- list(cif="strauss",par=c(beta=2,gamma=0,r=0.7),w=c(0,10,0,10))
   X3.strauss <- rmh(model=mod02,start=list(n.start=60),
                     control=list(nrep=nr,nverb=nv,iseed=c(42,17,69)))
   
   # Strauss process in a polygonal window.
   x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
   y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
   mod03 <- list(cif="strauss",par=c(beta=2000,gamma=0.6,r=0.07),
                w=owin(poly=list(x=x,y=y)))
   X4.strauss <- rmh(model=mod03,start=list(n.start=90),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss process in a polygonal window, conditioning on n = 80.
   X5.strauss <- rmh(model=mod03,start=list(n.start=90),
                     control=list(p=1,nrep=nr,nverb=nv))
   
   # Strauss process, starting off from X4.strauss, but with the
   # polygonal window replace by a rectangular one.  At the end,
   # the generated pattern is clipped to the original polygonal window.
   xxx <- X4.strauss
   xxx$window <- as.owin(c(0,1,0,1))
   X6.strauss <- rmh(model=mod03,start=list(x.start=xxx),
                     control=list(nrep=nr,nverb=nv))
   
   # Strauss with hardcore:
   mod04 <- list(cif="straush",par=c(beta=2,gamma=0.2,r=0.7,hc=0.3),
                w=c(0,10,0,10))
   X1.straush <- rmh(model=mod04,start=list(n.start=70),
                     control=list(nrep=nr,nverb=nv))
   
   # Another Strauss with hardcore (with a perhaps surprising result):
   mod05 <- list(cif="straush",par=c(beta=80,gamma=0.36,r=45,hc=2.5),
                w=c(0,250,0,250))
   X2.straush <- rmh(model=mod05,start=list(n.start=250),
                     control=list(nrep=nr,nverb=nv))
   
   # Pure hardcore (identical to X3.strauss).
   mod06 <- list(cif="straush",par=c(beta=2,gamma=1,r=1,hc=0.7),
                w=c(0,10,0,10))
   X3.straush <- rmh(model=mod06,start=list(n.start=60),
                     control=list(nrep=nr,nverb=nv,iseed=c(42,17,69)))
   
   # Soft core:
   par3 <- c(0.8,0.1,0.5)
   w    <- c(0,10,0,10)
   mod07 <- list(cif="sftcr",par=c(beta=0.8,sigma=0.1,kappa=0.5),
                w=c(0,10,0,10))
   X.sftcr <- rmh(model=mod07,start=list(n.start=70),
                  control=list(nrep=nr,nverb=nv))
   
   # Multitype Strauss:
   beta <- c(0.027,0.008)
   gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
   r    <- matrix(c(45,45,45,45),2,2)
   mod08 <- list(cif="straussm",par=list(beta=beta,gamma=gmma,radii=r),
                w=c(0,250,0,250))
   X1.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   
   # Multitype Strauss conditioning upon the total number
   # of points being 80:
   X2.straussm <- rmh(model=mod08,start=list(n.start=80),
                      control=list(p=1,ptypes=c(0.75,0.25),nrep=nr,
                                   nverb=nv))
   
   # Conditioning upon the number of points of type 1 being 60
   # and the number of points of type 2 being 20:
   X3.straussm <- rmh(model=mod08,start=list(n.start=c(60,20)),
                      control=list(fixall=TRUE,p=1,ptypes=c(0.75,0.25),
                                   nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore:
   rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
   mod09 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                iradii=r,hradii=rhc),w=c(0,250,0,250))
   X.straushm <- rmh(model=mod09,start=list(n.start=80),
                     control=list(ptypes=c(0.75,0.25),nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type:
   beta  <- c(0.27,0.08)
   tr3   <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
                         # log quadratic trend
   tr4   <- function(x,y){x <- x/250; y <- y/250;
                         exp(-0.6*x+0.5*y)}
                        # log linear trend
   mod10 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=c(0,250,0,250),
                 trend=list(tr3,tr4),tmax=list(1.5,1.65))
   X1.straushm.trend <- rmh(model=mod10,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),
                            nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type, given as images:
   x <- seq(0,250,length=51)
   xy <- expand.grid(x=x,y=x)
   i1 <- im(matrix(tr3(xy$x,xy$y),nrow=51),x,x)
   i2 <- im(matrix(tr4(xy$x,xy$y),nrow=51),x,x)
   mod11 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=c(0,250,0,250),
                 trend=list(i1,i2))
   X2.straushm.trend <- rmh(model=mod11,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),expand=1,
                            nrep=nr,nverb=nv))
   
   # Diggle, Gates, and Stibbard:
   mod12 <- list(cif="dgs",par=c(beta=3600,rho=0.08),w=c(0,1,0,1))
   X.dgs <- rmh(model=mod12,start=list(n.start=300),
                control=list(nrep=nr,nverb=nv))
   
   # Diggle-Gratton:
   mod13 <- list(cif="diggra",
                 par=c(beta=1800,kappa=3,delta=0.02,rho=0.04),
                 w=square(1))
   X.diggra <- rmh(model=mod13,start=list(n.start=300),
                   control=list(nrep=nr,nverb=nv))
   
   # Geyer:
   mod14 <- list(cif="geyer",par=c(beta=1.25,gamma=1.6,r=0.2,sat=4.5),
                 w=c(0,10,0,10))
   X1.geyer <- rmh(model=mod14,start=list(n.start=200),
                   control=list(nrep=nr,nverb=nv))
   
   # Geyer; same as a Strauss process with parameters
   # (beta=2.25,gamma=0.16,r=0.7):
   
   mod15 <- list(cif="geyer",par=c(beta=2.25,gamma=0.4,r=0.7,sat=10000),
                 w=c(0,10,0,10))
   X2.geyer <- rmh(model=mod15,start=list(n.start=200),
                   control=list(nrep=nr,nverb=nv))
   
   mod16 <- list(cif="geyer",par=c(beta=8.1,gamma=2.2,r=0.08,sat=3))
   data(redwood)
   X3.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(periodic=TRUE,nrep=nr,nverb=nv))
   
   # Geyer, starting from the redwood data set, simulating
   # on a torus, and conditioning on n:
   X4.geyer <- rmh(model=mod16,start=list(x.start=redwood),
                   control=list(p=1,periodic=TRUE,nrep=nr,nverb=nv))
   
   # Strauss with trend
   tr <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
   beta <- 0.3
   gmma <- 0.5
   r    <- 45
   mod17 <- list(cif="strauss",par=c(beta=beta,gamma=gmma,r=r),w=c(0,250,0,250),
                 trend=tr3,tmax=1.5)
   X1.strauss.trend <- rmh(model=mod17,start=list(n.start=90),
                           control=list(nrep=nr,nverb=nv))
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
