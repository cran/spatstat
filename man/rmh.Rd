\name{rmh}
\alias{rmh}
\title{Simulate point patterns using the Metropolis-Hastings algorithm.}
\description{
Simulates point patterns from a (still small) range of point process models.
}
\usage{
rmh(cif,par,w=NULL,ntypes=1,ptypes=NULL,tpar=NULL,n.start=NULL,
        x.start=NULL,fixall=FALSE,expand=NULL,periodic=F,nrep=1e6,
        p=0.9,q=0.5,iseed=NULL,nverb=0)
}
\arguments{
  \item{cif}{A character string giving the \bold{name} of a
    Fortran subroutine which will calculate the value of the
    conditional intensity function for the model.  So far the only
    permissible values are \code{strauss} (Strauss process),
    \code{straush} (Strauss process with hardcore), \code{sftcr}
    (Softcore process), \code{straussm} (multitype Strauss process),
    \code{straushm} (multitype Strauss process with hardcore),
    \code{dig1}, \code{dig2} (see the 2nd of the \bold{References}),
    and \code{geyer} (see the 3rd of the \bold{References}).

    It is intended that more options will be added in the future.
    The very brave user \bold{could} try to add her own.  Note that
    in addition to writing Fortran or ratfor code for the new
    conditional intensity function, the user would have to modify the
    code in the files \code{cif.r} and \code{rmh.R} appropriately.
    (And of course re-install the \code{spatstat} package so as to
    update the dynamically loadable shared object \code{spatstat.so}.)
  }

  \item{par}{A set of parameter values appropriate to the conditional
    intensity function being invoked.  (See \bold{Details}.)
  }

  \item{w}{A specification of a window in which the pattern is
    to be generated.  If specified, it must be in a form which can be
    coerced to an object of class \code{owin} by
    \code{\link{as.owin}()}.  Note that for non-rectangular windows
    we cannot (as yet) generated a pattern from a process which
    exists only within that window.  The domain of the theoretical
    process must (at present) be at least as large as the ``enclosing
    box'' of the window.
  }

  \item{ntypes}{The number of distinct types to be used if
    the simulated pattern is to be a multitype (marked) point
    pattern.
  }
  \item{ptypes}{A vector of probabilities (of length \code{ntypes}
    and summing to 1) to be used in assigning a random type to a new
    point.  Defaults to a vector each of whose entries is
    1/\code{ntypes}.  Convergence of the simulation algorithm should
    be improved if \code{ptypes} is close to the relative frequencies
    of the types which will result from the simulation.
  }
  \item{tpar}{A vector (or possibly a list if the pattern is multitype)
    specifying the coefficients of a polynomial for a log polynomial
    trend.  The coefficients (corresponding to each type, for a
    multitype process) must be given in the order of the terms
    \eqn{x, y, x^2, xy, y^2, x^3, x^2y, xy^2, y^3, \ldots, y^n} where
    \eqn{n} is the degree of the trend.  (Thus the possible lengths
    of a sequence of coefficients are 0, 2, 5, 9, 14, \eqn{\ldots}).

    If the process is unmarked \code{tpar} must be a vector
    consisting solely of the coefficients of the trend polynomial in
    the correct order.  The degree of the polynomial will be
    calculated from the length of \code{tpar}.

    If the process is multitype and if \code{tpar} is given as a vector,
    then \code{ntypes} and and the \bold{degrees} of each of the
    polynomial should comprise the first 1+\code{ntypes} entries
    of \code{tpar}.  After that, the coefficients corresponding to
    each mark should be catenated successively.  Note that some
    of the degrees are allowed to be 0.

    If \code{tpar} is given as a list, the \eqn{i^{th}}{i-th}
    component of this list should consist of the coefficients for the
    trend corresponding to the \eqn{i^{th}}{i-th} type.  If there is
    to be no trend for the \eqn{i^{th}}{i-th} type, then the
    \eqn{i^{th}}{i-th} component of the list should be NULL (but must
    be present in the list as the NULL object.)
  }
 
  \item{n.start}{The number of ``initial'' points to be randomly
    (uniformly) generated in the window \code{w}, or a vector (of
    length \code{ntypes}) giving the number of points of each type to
    be generated.  These uniformly generated points form an initial
    state or configuration for the Metropolis-Hastings algorithm.
    Precisely \bold{one} of \code{n.start} or \code{x.start} must
    be specified.

    A vector-valued \code{n.start} is meaningful \bold{only}
    if \code{p} (the probability of a shift as opposed to a birth or
    death) is equal to 1 (so that we are conditioning on the number
    of points).  When \eqn{p < 1}, if \code{n.start} is vector valued
    then it is effectively replaced by its sum.

    The resulting set of uniformly generated points gives the
    Metropolis-Hastings algorithm an initial state from which to
    start.  (Actually, when \eqn{p < 1}, the number \code{n.start}
    gets multiplied by the ratio of the area of the enclosing box for
    the window to the area of the window, and then by the factor
    \code{expand}.  Then that many points are uniformly generated in
    the \bold{expanded} window; see below.)  The value of
    \code{n.start} should be roughly equal to (an educated guess at)
    the expected number of points which will be generated inside the
    window.
   }

  \item{x.start}{A point pattern (an object of class \code{ppp}),
    which provides an alternative means of specifying the initial
    ``state'' or configuration for the Metropolis-Hastings algorithm.
    The simulated pattern is constructed in the enclosing rectangle
    of \code{x.start$w}.  No expansion takes place.  (The argument
    \code{expand} is forced to equal 1.)  At the end of the
    simulation, the resulting pattern is clipped to the window
    \code{w} if this is given, otherwise to \code{x.start$w}.
  }
  
  \item{fixall}{A logical scalar specifying whether to condition on
   the number of points of each type.  Meaningful only if a marked
   process is being simulated, and if \eqn{p = 1}.  A warning message
   is given if \code{fixall} is set equal to \code{TRUE} when it is
   not meaningful.
  }
 
  \item{expand}{The factor by which the enclosing box of the window
    \code{w} is to be expanded in order to better approximate the
    simulation of a process existing in the whole plane, rather than
    just in the enclosing box.  If \code{expand} equals 1, then we
    are simulating the latter (unless \code{periodic} in \code{TRUE};
    see below).  The larger \code{expand} is, the better we
    approximate the former.  Note that any value of \code{expand}
    smaller than 1 is treated as if it were 1.

    The area of the expanded window is equal to \code{expand} times
    the area of the enclosing box; width and height are stretched
    proportionately.  Points are generated by the Metropolis-Hastings
    algorithm in the expanded window, and then ``clipped'' down to
    the original window \code{w} (or to \code{x.start$w}) when the
    algorithm has finished.  The argument \code{expand} defaults to 2
    if \code{periodic} is \code{FALSE} and \eqn{p < 1} and to 1 if
    \code{periodic} is \code{TRUE} or if \eqn{p < 1}, or if the
    starting configuration is specified via \code{x.start}.  Trying to
    set \code{expand} greater than 1 when \code{periodic} is
    \code{TRUE} or \eqn{p = 1} generates an error.  A specified value
    of \code{expand} is simply ignored if \code{x.start} is given.
  }
  \item{periodic}{A logical scalar; if \code{periodic} is \code{TRUE}
    we simulate a process on the torus formed by identifying opposite
    edges of the (rectangular) window.  If \code{periodic} is \code{TRUE}
    and the window (\code{x.start$w} if \code{x.start} is specified;
    otherwise \code{w}) is not rectangular, an error is given.
  }

  \item{nrep}{The number of repetitions or steps (changes of state)
    to be made by the Metropolis-Hastings algorithm.  It should
    be large.
  }

  \item{p}{The probability of of proposing a ``shift'' (as opposed
    to a birth or death) in the Metropolis-Hastings algorithm.  If
    \eqn{p = 1} then we do nothing but shifts, whence the number of
    points never changes, i.e. we are simulating conditionally upon
    the number of points.  In this case, for multitype processes, we
    also have the option of simulating conditionally upon the number
    of points of each type; this is effected by specifying
    \code{fixall=TRUE}.  In this case, \code{n.start} must be a
    vector, of length \code{ntypes}, whose entries are these
    numbers.

    We can only condition on the number of points if the simulation
    takes place in the original window (as opposed to taking place in
    a larger window and then being clipped to the original).  Hence,
    if \eqn{p = 1} then \code{expand} defaults to 1 and it is an
    error to specify a value of \code{expand} which is greater than
    1.  Likewise, if \eqn{p = 1}, it is (currently) an error to specify
    anything other than a rectangular window.
  }

  \item{q}{The probability of proposing a death (rather than a birth)
    given that birth/death has been chosen over shift.  This is of
    course ignored if \code{p} is equal to 1.
  }

  \item{iseed}{A vector equal to a triple of integers to be used
    as seeds to the random number generating procedure.
    If unspecified these are themselves generated, on the interval
    from 1 to 1 million, using the function \code{sample()}.
  }

  \item{nverb}{An integer specifying how often ``progress reports''
    (which consist simply of the number of repetitions completed)
    should be printed out.  If nverb is left at 0, the default,
    the simulation proceeds silently.
  }
}

\value{ A list of class \code{ppp} with the usual components

  \item{window}{The window \code{w} in the form of
    an object of class \code{owin}.
  }

  \item{n}{The number of generated points (in the window \code{w}).}

  \item{x}{The x-coordinates of the generated points.}

  \item{y}{The y-coordinates of the generated points.}

  \item{marks}{The marks of the generated points.  Remember that this
   component is a \bold{factor}.  (This component is present only
   when the model is a multitype point process.)}

   In addition the list returned has a component \code{info}
   consisting of arguments supplied to the function (or default
   values of arguments which were not explicitly supplied).  These
   are given so that it is possible to reconstruct exactly the manner
   in which the pattern was generated.  The components of \code{info}
   are: \code{cif}, \code{par}, \code{tpar}, \code{n.start},
   \code{x.start.name}, \code{nrep}, \code{p}, \code{q},
   \code{expand}, \code{periodic}, and \code{iseed}.

   Note that if \code{x.start} is specified only its \bold{name} is
   preserved in the list \code{info}.

}

\details{
  The parameter vector (or list) \code{par} should be
  as follows, for each of the available conditional intensity
  functions:
  \describe{
    \item{strauss:}{
      (Strauss process.) A vector with components
      \eqn{(\beta,\gamma,r)}{beta,gamma,r} which are respectively the
      ``base'' intensity, the pair-wise interaction parameter and the
      interaction radius.  Note that \eqn{\gamma}{gamma} must be less
      than or equal to 1.
    }
    \item{straush:}{
      (Strauss process with hardcore.) A vector with entries
      \eqn{(\beta,\gamma,r,r_{hc})}{beta,gamma,r,r_hc} where
      \eqn{\beta}{beta}, \eqn{\gamma}{gamma}, and \eqn{r} are as for
      the Strauss process, and \eqn{r_{hc}}{r_hc} is the hardcore
      radius.  Of course \eqn{r_{hc}}{r_hc} must be less than
      \eqn{r}.
    }
    \item{sftcr:}{
      (Softcore process.) A vector with components
      \eqn{(\beta,\sigma,\kappa)}{beta,sigma,kappa}.  Again
      \eqn{\beta}{beta} is a ``base'' intensity.  The pairwise
      interaction between two points \eqn{u \neq v}{u != v} is
      \deqn{
	\exp \left \{
	- \left ( \frac{\sigma}{||u-v||} \right )^{2/\kappa}
	\right \}
      }{-(sigma/||u-v||)^(2/kappa)}
      Note that it is necessary that \eqn{0 < \kappa < 1}{0 < kappa <1}.
    }
    \item{straussm:}{
      (Multitype Strauss process.) Here \code{par} is
      best given as a list with components
      \describe{
	\item{beta:}{A vector of ``base'' intensities, one for each possible
	  type.}
	\item{gamma:}{A \bold{symmetric} matrix of interaction parameters,
	  with \eqn{\gamma_{ij}}{gamma_ij} pertaining to the interaction between
	  type \eqn{i} and type \eqn{j}.
	}
	\item{r:}{A \bold{symmetric} matrix of interaction radii, with
	  \eqn{r_{ij}}{r_ij} pertaining to the interaction between type
	  \eqn{i} and type \eqn{j}.}
      }
      If \code{par} is to be given as a vector, then this vector should
      consist of the vector \code{beta}, catenated with the upper
      triangle of the matrix \code{gamma} strung out in \bold{row}
      order, and finally catenated with the upper triangle of the
      matrix \code{r} likewise strung out in \bold{row} order.
    }
    \item{straushm:}{
      (Multitype Strauss process with hardcore.) The
      parameters \code{par} are much as for \code{straussm} except
      that there is an extra component \code{rhc} which is the matrix
      of hardcore radii.  Note that in the vector form, the hardcore
      radii must be catenated \bold{after} the interaction radii.
    }
    \item{dig1:}{
      (See the 2nd of the \bold{References}.)
      Process with pairwise interaction function
      \deqn{e(t) = \sin^2(\pi t/2\rho)}{e(t) = sin^2((pi t)/(2 \rho))}
      for \eqn{t < \rho}{t < rho}, and equal to 1 for \eqn{t \geq
	\rho}{t >= rho}.  The parameters are \eqn{\beta}{beta} and
      \eqn{\rho}{rho}.
    }
    \item{dig2:}{
      (See the 2nd of the \bold{References}.)
      Process with pairwise interaction function \eqn{e(t)} equal to 0
      for \eqn{t < \delta}{t < delta}, equal to
      \deqn{((t-\delta)/(\rho-\delta))^{\kappa}}{
	((t-delta)/(rho-delta))^kappa} for \eqn{\delta \leq t <
	\rho}{delta <= t < rho}, and equal to 1 for \eqn{\delta \geq
	\rho}{delta >= rho}.  Note that here we use the symbol
      \eqn{\kappa}{kappa} where Diggle, Gates, and Stibbard use
      \eqn{\beta}{beta} since we reserve the symbol \eqn{\beta}{beta}
      for an intensity parameter.
      
      The complete parameter set is \eqn{\beta}{beta}, \eqn{\kappa}{kappa},
      \eqn{\delta}{delta} and \eqn{\rho}{rho}.
    }
    \item{geyer}{
      (See the 3rd of the \bold{References}.) Geyer's ``saturation''
      point process model.  This model is ``like a Strauss model,
      but with an upper bound to the number of \eqn{r}-close neighbors
      of any point.''

      More explicitly, a saturation point process with interaction
      radius \eqn{r}, saturation threshold \eqn{s}, and  parameters
      \eqn{\beta}{beta} and \eqn{\gamma}{gamma}, is the point process
      in which each point \eqn{x_i} in the pattern \eqn{X} contributes
      a factor
      \deqn{\beta \gamma^{\min(s, t(x_i,X))}}{beta gamma^min(s,t(x_i,X))}
      to the probability density of the point pattern, where
      \eqn{t(x_i,X)} denotes the number of ``\eqn{r}-close neighbours''
      of \eqn{x_i} in the pattern \eqn{X}.
    }
  }
}

\references{
   Baddeley, Adrian, and Turner, Rolf.  ``Practical maximum
   pseudolikelihood for spatial point patterns.'' Australian and New
   Zealand Journal of Statistics, vol. 42, 2000, pp. 283 -- 322.

   Diggle, Peter J., Gates, David J., and Stibbard, Alyson.  ``A
   nonparametric estimator for pairwise-interaction point processes.''
   Biometrika, vol. 74, 1987, pp. 763 -- 770.

   Geyer, C.J. (1999) Likelihood Inference for Spatial Point
   Processes. Chapter 3 in  O.E. Barndorff-Nielsen, W.S. Kendall and
   M.N.M. Van Lieshout (eds) Stochastic Geometry: Likelihood and
   Computation, Chapman and Hall / CRC,  Monographs on Statistics and
   Applied Probability, number 80. Pages 79-140
}

\note{It is possible to simulate conditionally upon the number of
points, or in the case of multitype processes, upon the number of
points of each type.  To condition upon the total number of points,
set \code{p} (the probability of a shift) equal to 1, and specify
\code{n.start} to be a scalar (as usual).  To condition upon the
number of points of each type, set \code{p} equal to 1, \code{fixall}
equal to \code{TRUE}, and specify \code{n.start} to be a vector (of
length \code{ntypes}).

In these circumstances
	\itemize{
	\item The value of \code{expand} must be equal to 1; it
	defaults to 1, and it is an error to specify a value larger
	than 1.
	\item The window in which the pattern is being simulated
	must (currently) be rectangular.
	\item The resulting simulated pattern will have precisely
	the number of points (of each type) specified by \code{n.start}.
	}
}

\section{Warnings}{
There is never a guarantee that the Metropolis-Hastings algorithm
has converged to the steady state.

If \code{x.start} is specified then \code{expand} is set equal
to 1 and simulation takes place in the enclosing rectangle of
\code{x.start$w}.  Any specified value for \code{expand} is
simply ignored.

If trends are specified, make sure that the lengths of the
vectors of coefficients in \code{tpar} make sense.  For multitype
processes, it is probably safest to specify \code{tpar} as a
list.  But make sure that, even if there is to be no trend
corresponding to a particular type, there is still a component
(a NULL component) for that type, in the list.

Note that if \code{tpar} is given as an atomic vector, for
multitype processes, no checking is (or, realistically, can be)
done to make sure that the degrees supplied are sensible.
}
\seealso{
  \code{\link{ppp}},
  \code{\link{mpl}},
  \code{\link{Strauss}},
  \code{\link{Softcore}},
  \code{\link{StraussHard}},
  \code{\link{MultiStrauss}},
  \code{\link{MultiStraussHard}}
}
\examples{
   \dontrun{library(spatstat)}
   \dontrun{
    par11 <- c(2,0.2,0.7)
    w     <- c(0,10,0,10)
    X1.strauss <- rmh("strauss",par=par11,w=w,n.start=80,
                      nrep=1e5,nverb=5000)

    # Conditioning on n = 80:
    X2.strauss <- rmh("strauss",par=par11,w=w,n.start=80,
                      p=1,nrep=1e5,nverb=5000)
    
    par12 <- c(2000,0.6,0.07)
    x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
    y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
    w     <- owin(poly=list(x=x,y=y))
    X3.strauss <- rmh("strauss",par=par12,w=w,n.start=90,
                      nrep=1e5,nverb=5000)
    
    # Pure hardcore:
    par13 <- c(2,0,0.7)
    w     <- c(0,10,0,10)
    X4.strauss <- rmh("strauss",par=par13,w=w,n.start=60,
                      nrep=1e5,nverb=5000,iseed=c(42,17,69))

    # Another Strauss process, starting off from X3.strauss (but
    # with a rectangular window):
    xxx <- X3.strauss
    xxx$w <- as.owin(c(0,1,0,1))
    X5.strauss <- rmh("strauss",par=par12,w=w,x.start=xxx,
                      nrep=1e5,nverb=5000)
    
    par21 <- c(2,0.2,0.7,0.3)
    w     <- c(0,10,0,10)
    X1.straush <- rmh("straush",par=par21,w=w,n.start=70,
                      nrep=1e5,nverb=5000)
    
    par22 <- c(80,0.36,45,2.5)
    w     <- c(0,250,0,250)
    X2.straush <- rmh("straush",par=par22,w=w,n.start=160,
                      nrep=1e5,nverb=5000)
    
    # Pure hardcore (identical to X3.strauss).
    par23 <- c(2,1,1,0.7)
    w     <- c(0,10,0,10)
    X3.straush <- rmh("straush",par=par23,w=w,n.start=60,
                      nrep=1e5,nverb=5000,iseed=c(42,17,69))
    
    par3 <- c(0.8,0.1,0.5)
    w    <- c(0,10,0,10)
    X.sftcr <- rmh("sftcr",par=par3,w=w,n.start=70,nrep=1e5,
                   nverb=5000)
    
    beta <- c(0.027,0.008)
    gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
    r    <- matrix(c(45,45,45,45),2,2)
    par4 <- list(beta=beta,gamma=gmma,r=r)
    w    <- c(0,250,0,250)
    pm   <- c(0.75,0.25)
    X1.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,
                      ptypes=pm,n.start=80,nrep=1e5,nverb=5000)
    
    # Conditioning upon the total number of points being 80:
    X2.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,
                      ptypes=pm,n.start=80,p=1,nrep=1e5,nverb=5000)
    
    # Conditioning upon the number of points of type 1 being 60
    # and the number of points of type 2 being 20:
    X3.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,n.start=c(60,20),
                       fixall=TRUE,nrep=1e5,p=1,nverb=5000)

    rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
    par5 <- list(beta=beta,gamma=gmma,r=r,rhc=rhc)
    X.straushm <- rmh("straushm",par=par5,w=w,ntypes=2,
                      ptypes=pm,n.start=80,nrep=1e5,nverb=5000)
    
    beta  <- c(0.0027,0.08)
    par6  <- list(beta=beta,gamma=gmma,r=r,rhc=rhc)
    tpar1 <- c(0.02,0.004,-0.0004,0.004,-0.0004) # Coefs. for log quadratic
    tpar2 <- c(-0.06,0.05)                       # and log linear trends.
    w     <- c(0,250,0,250)
    pm    <- c(0.75,0.25)
    X1.straushm.trend <- rmh("straushm",par=par6,w=w,ntypes=2,
                             ptypes=pm,tpar=list(tpar1,tpar2),n.start=350,
                             nrep=1e5,nverb=5000,iseed=c(42,17,69))
    
    # Identical to X1.straushm.trend; tpar given in vector form.
    tpar <- c(2,2,1,tpar1,tpar2)
    X2.straushm.trend <- rmh("straushm",par=par6,w=w,ntypes=2,
                             ptypes=pm,tpar=tpar,n.start=350,
                             nrep=1e5,nverb=5000,iseed=c(42,17,69))
    par7 <- c(3600,0.08)
    w    <- c(0,1,0,1)
    X.dig1 <- rmh("dig1",par=par7,w=w,n.start=300,nrep=1e5,nverb=5000)
    
    par8 <- c(1800,3,0.02,0.04)
    X.dig2 <- rmh("dig2",par=par8,w=w,n.start=300,nrep=1e5,nverb=5000)
    
    par9 <- c(1.25,1.6,0.2,4.5)
    w    <- c(0,10,0,10)
    X1.geyer <- rmh("geyer",par=par9,w=w,n.start=200,nrep=1e5,nverb=5000)

    # Same as a Strauss process with parameters (2.25,0.16,0.7).
    par10 <- c(2.25,0.4,0.7,10000)
    X2.geyer <- rmh("geyer",par=par10,w=w,n.start=70,nrep=1e5,nverb=5000)

    # Starting from the redwood data set, and simulating on a torus:
    par11 <- c(8.1,2.2,0.08,3)
    data(redwood)
    X3.geyer <- rmh("geyer",par=par11,x.start=redwood,periodic=TRUE,
                    nrep=1e5,nverb=5000)

    # As above, conditioning on n:
    X4.geyer <- rmh("geyer",par=par11,x.start=redwood,periodic=TRUE,
                    p=1,nrep=1e5,nverb=5000)
   }
   \testonly{
    par11 <- c(2,0.2,0.7)
    w     <- c(0,10,0,10)
    X1.strauss <- rmh("strauss",par=par11,w=w,n.start=80,
                      nrep=100,nverb=50)

    # Conditioning on n = 80:
    X2.strauss <- rmh("strauss",par=par11,w=w,n.start=80,
                      p=1,nrep=100,nverb=50)
    
    par12 <- c(2000,0.6,0.07)
    x     <- c(0.55,0.68,0.75,0.58,0.39,0.37,0.19,0.26,0.42)
    y     <- c(0.20,0.27,0.68,0.99,0.80,0.61,0.45,0.28,0.33)
    w     <- owin(poly=list(x=x,y=y))
    X3.strauss <- rmh("strauss",par=par12,w=w,n.start=90,
                      nrep=100,nverb=50)
    
    # Pure hardcore:
    par13 <- c(2,0,0.7)
    w     <- c(0,10,0,10)
    X4.strauss <- rmh("strauss",par=par13,w=w,n.start=60,
                      nrep=100,nverb=50,iseed=c(42,17,69))

    # Another Strauss process, starting off from X3.strauss (but
    # with a rectangular window):
    xxx <- X3.strauss
    xxx$w <- as.owin(c(0,1,0,1))
    X5.strauss <- rmh("strauss",par=par12,w=w,x.start=xxx,
                      nrep=100,nverb=50)
    
    par21 <- c(2,0.2,0.7,0.3)
    w     <- c(0,10,0,10)
    X1.straush <- rmh("straush",par=par21,w=w,n.start=70,
                      nrep=100,nverb=50)
    
    par22 <- c(80,0.36,45,2.5)
    w     <- c(0,250,0,250)
    X2.straush <- rmh("straush",par=par22,w=w,n.start=160,
                      nrep=100,nverb=50)
    
    # Pure hardcore (identical to X3.strauss).
    par23 <- c(2,1,1,0.7)
    w     <- c(0,10,0,10)
    X3.straush <- rmh("straush",par=par23,w=w,n.start=60,
                      nrep=100,nverb=50,iseed=c(42,17,69))
    
    par3 <- c(0.8,0.1,0.5)
    w    <- c(0,10,0,10)
    X.sftcr <- rmh("sftcr",par=par3,w=w,n.start=70,nrep=100,
                   nverb=50)
    
    beta <- c(0.027,0.008)
    gmma <- matrix(c(0.43,0.98,0.98,0.36),2,2)
    r    <- matrix(c(45,45,45,45),2,2)
    par4 <- list(beta=beta,gamma=gmma,r=r)
    w    <- c(0,250,0,250)
    pm   <- c(0.75,0.25)
    X1.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,
                      ptypes=pm,n.start=80,nrep=100,nverb=50)
    
    # Conditioning upon the total number of points being 80:
    X2.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,
                      ptypes=pm,n.start=80,p=1,nrep=100,nverb=50)
    
    # Conditioning upon the number of points of type 1 being 60
    # and the number of points of type 2 being 20:
    X3.straussm <- rmh("straussm",par=par4,w=w,ntypes=2,n.start=c(60,20),
                       fixall=TRUE,nrep=100,p=1,nverb=50)
    
    rhc  <- matrix(c(9.1,5.0,5.0,2.5),2,2)
    par5 <- list(beta=beta,gamma=gmma,r=r,rhc=rhc)
    X.straushm <- rmh("straushm",par=par5,w=w,ntypes=2,
                      ptypes=pm,n.start=80,nrep=100,nverb=50)
    
    beta  <- c(0.0027,0.08)
    par6  <- list(beta=beta,gamma=gmma,r=r,rhc=rhc)
    tpar1 <- c(0.02,0.004,-0.0004,0.004,-0.0004) # Coefs. for log quadratic
    tpar2 <- c(-0.06,0.05)                       # and log linear trends.
    w     <- c(0,250,0,250)
    pm    <- c(0.75,0.25)
    X1.straushm.trend <- rmh("straushm",par=par6,w=w,ntypes=2,
                             ptypes=pm,tpar=list(tpar1,tpar2),n.start=350,
                             nrep=100,nverb=50,iseed=c(42,17,69))
    
    # Identical to X1.straushm.trend; tpar given in vector form.
    tpar <- c(2,2,1,tpar1,tpar2)
    X2.straushm.trend <- rmh("straushm",par=par6,w=w,ntypes=2,
                             ptypes=pm,tpar=tpar,n.start=350,
                             nrep=100,nverb=50,iseed=c(42,17,69))
    par7 <- c(3600,0.08)
    w    <- c(0,1,0,1)
    X.dig1 <- rmh("dig1",par=par7,w=w,n.start=300,nrep=100,nverb=50)
    
    par8 <- c(1800,3,0.02,0.04)
    X.dig2 <- rmh("dig2",par=par8,w=w,n.start=300,nrep=100,nverb=50)
    
    par9 <- c(1.25,1.6,0.2,4.5)
    w    <- c(0,10,0,10)
    X1.geyer <- rmh("geyer",par=par9,w=w,n.start=200,nrep=100,nverb=50)

    # Same as a Strauss process with parameters (2.25,0.16,0.7).
    par10 <- c(2.25,0.4,0.7,10000)
    X2.geyer <- rmh("geyer",par=par10,w=w,n.start=70,nrep=100,nverb=50)

    # Starting from the redwood data set, and simulating on a torus:
    par11 <- c(8.1,2.2,0.08,3)
    data(redwood)
    X3.geyer <- rmh("geyer",par=par11,x.start=redwood,periodic=TRUE,
                    nrep=100,nverb=50)

    # As above, conditioning on n:
    X4.geyer <- rmh("geyer",par=par11,x.start=redwood,periodic=TRUE,
                    p=1,nrep=100,nverb=50)
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
