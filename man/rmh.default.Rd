\name{rmh.default}
\alias{rmh.default}
\title{Simulate Point Process Models using the Metropolis-Hastings Algorithm.}
\description{
  Generates a random point pattern, simulated from
  a chosen point process model, using the Metropolis-Hastings
  algorithm. 
}

\usage{\method{rmh}{default}(model, start, control, verbose=TRUE, \dots) }

\arguments{
  \item{model}{Data specifying the point process model
    that is to be simulated.
  }
  \item{start}{Data determining the initial state of
    the algorithm.
  }
  \item{control}{Data controlling the iterative behaviour
    and termination of the algorithm.
  }
  \item{verbose}{Logical flag indicating whether to print progress
    reports.
  }
  \item{\dots}{Further arguments which will be passed to
    any trend functions in \code{model}.
  }
}

\value{A point pattern (an object of class \code{"ppp"}, see
  \code{\link{ppp.object}}).

  The returned value has an attribute \code{info} containing
  modified versions of the arguments
  \code{model}, \code{start}, and \code{control} which together specify
  the exact simulation procedure.
}

\details{
  This function generates simulated realisations from any of a range of
  spatial point processes, using the Metropolis-Hastings algorithm.
  It is the default method for the generic function \code{\link{rmh}}.

  This function executes a Metropolis-Hastings algorithm
  with birth, death and shift proposals as described in
  Geyer and Moller (1994).

  The argument \code{model} specifies the point process model to be
  simulated. It is either a list, or an object of class
  \code{"rmhmodel"}, with the following components:

  \describe{
    \item{cif}{A character string specifying the choice of
      interpoint interaction for the point process.
    }
    \item{par}{
      Parameter values for the conditional
      intensity function.
    }
    \item{w}{
      (Optional) window in which the pattern is
      to be generated. An object of class \code{"owin"},
      or data acceptable to \code{\link{as.owin}}.
    }
    \item{trend}{
      Data specifying the spatial trend in the model, if it has a trend.
      This may be a function, a pixel image (of class \code{"im"}),
      (or a list of functions or images if the model
      is multitype).
      
      If the trend is a function or functions,
      any auxiliary arguments \code{...} to \code{rmh.default}
      will be passed to these functions, which
      should be of the form \code{function(x, y, ...)}.
    }
    \item{types}{
      List of possible types, for a multitype point process.
    }
  }
  For full details of these parameters, see \code{\link{rmhmodel}}.

  
  The argument \code{start} determines the initial state of the
  Metropolis-Hastings algorithm. It is either a list, or an object
  of class \code{"rmhstart"}, with the following components:

  \describe{
    \item{n.start}{
      Number of points in the initial point pattern.
      A single integer, or a vector of integers giving the
      numbers of points of each type in a multitype point pattern.
      Incompatible with \code{x.start}.
    }
    \item{x.start}{
      Initial point pattern configuration.
      Incompatible with \code{n.start}.

      \code{x.start} may be a point pattern (an
      object of class \code{"ppp"}), or data which can be coerced
      to this class by \code{\link{as.ppp}},  or an object with
      components \code{x} and \code{y}, or a two-column matrix.
      In the last two cases, the window for the pattern is determined
      by \code{model$w}.
      In the first two cases, if \code{model$w} is also present,
      then the final simulated pattern will be clipped to
      the window \code{model$w}.
    }
    \item{iseed}{
      (Optional) seed for random number generator.
      A triple of integers. Should not be specified
      in normal usage!
    }
  }
  For full details of these parameters, see \code{\link{rmhstart}}.

  The third argument \code{control} controls the simulation
  procedure, iterative behaviour, and termination of the
  Metropolis-Hastings algorithm. It is either a list, or an
  object of class \code{"rmhcontrol"}, with components:
  \describe{
    \item{p}{The probability of proposing a ``shift''
      (as opposed to a birth or death) in the Metropolis-Hastings
      algorithm.
    }
    \item{q}{The conditional probability of proposing a death
      (rather than a birth)
      given that birth/death has been chosen over shift.  
    }
    \item{nrep}{The number of repetitions or iterations
      to be made by the Metropolis-Hastings algorithm.  It should
      be large.
    }
    \item{expand}{
      Either a numerical expansion factor, or
      a window (object of class \code{"owin"}). Indicates that
      the process is to be simulated on a larger domain than the
      original data window \code{w}, then clipped to \code{w}
      when the algorithm has finished.

      If the model has a trend, then in order for expansion to
      be feasible, the trend must be given either as a function,
      or an image whose bounding box is large enough to contain
      the expanded window.
    }
    \item{periodic}{A logical scalar; if \code{periodic} is \code{TRUE}
      we simulate a process on the torus formed by identifying
      opposite edges of a rectangular window.  
    }
    \item{ptypes}{A vector of probabilities (summing to 1) to be used
      in assigning a random type to a new point.
    }
    \item{fixall}{A logical scalar specifying whether to condition on
      the number of points of each type.
    }
    \item{nverb}{An integer specifying how often ``progress reports''
      (which consist simply of the number of repetitions completed)
      should be printed out.  If nverb is left at 0, the default,
      the simulation proceeds silently.
    }
  }
  For full details of these parameters, see \code{\link{rmhcontrol}}.

  It is possible to simulate the model conditionally upon the number of
  points, or in the case of multitype processes, upon the number of
  points of each type.  To condition upon the total number of points,
  set \code{control$p} (the probability of a shift) equal to 1.
  The number of points is then determined by the starting state, which
  may be specified either by setting \code{start$n.start} to be a
  scalar, or by setting the initial pattern \code{start$x.start}.

  To condition upon the
  number of points of each type, set \code{control$p} equal to 1
  and \code{control$fixall} to be \code{TRUE}.
  The number of points is then determined by the starting state, which
  may be specified either by setting \code{start$n.start} to be an
  integer vector, or by setting the initial pattern \code{start$x.start}.
  
  When we simulate conditionally, no expansion of the window is possible.
}
     
\references{
   Baddeley, A. and Turner, R. (2000) Practical maximum
   pseudolikelihood for spatial point patterns.
   \emph{Australian and New Zealand Journal of Statistics}
   \bold{42}, 283 -- 322.

   Diggle, P. J. (2003) \emph{Statistical Analysis of Spatial Point
   Patterns} (2nd ed.) Arnold, London.

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

\section{Warnings}{

There is never a guarantee that the Metropolis-Hastings algorithm
has converged to its limiting distribution.

If \code{start$x.start} is specified then \code{expand} is set equal to 1
and simulation takes place in \code{x.start$window}.  Any specified
value for \code{expand} is simply ignored.

The presence of both a component \code{w} of \code{model} and a
non-null value for \code{x.start$window} makes sense ONLY if \code{w}
is contained in \code{x.start$window}.  

For multitype processes make sure that, even if there is to be no
trend corresponding to a particular type, there is still a component
(a NULL component) for that type, in the list.
}

\seealso{
  \code{\link{rmh}},
  \code{\link{rmh.ppm}},
  \code{\link{rStrauss}},
  \code{\link{ppp}},
  \code{\link{ppm}},
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
  the \pkg{spatstat} package so as to update the dynamically
  loadable shared object \code{spatstat.so}.)

  Note that the \code{lookup} conditional intensity function
  permits the simulation (in theory, to any desired degree
  of approximation) of any pairwise interaction process for
  which the interaction depends only on the distance between
  the pair of points.
}
\examples{
   \dontrun{
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
   X3.strauss <- rmh(model=mod02,start=list(n.start=60,iseed=c(42,17,69)),
                     control=list(nrep=nr,nverb=nv))
   
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
   X3.straush <- rmh(model=mod06,start=list(n.start=60, iseed=c(42,17,69)),
                     control=list(nrep=nr,nverb=nv))
   
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
                 trend=list(tr3,tr4))
   X1.straushm.trend <- rmh(model=mod10,start=list(n.start=350),
                            control=list(ptypes=c(0.75,0.25),
                            nrep=nr,nverb=nv))
   
   # Multitype Strauss hardcore with trends for each type, given as images:
   bigwin <- square(250)
   i1 <- as.im(tr3, bigwin)
   i2 <- as.im(tr4, bigwin)   
   mod11 <- list(cif="straushm",par=list(beta=beta,gamma=gmma,
                 iradii=r,hradii=rhc),w=bigwin,
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

   # Lookup (interaction function h_2 from page 76, Diggle (2003)):
      r <- seq(from=0,to=0.2,length=101)[-1] # Drop 0.
      h <- 20*(r-0.05)
      h[r<0.05] <- 0
      h[r>0.10] <- 1
      mod17 <- list(cif="lookup",par=list(beta=4000,h=h,r=r),w=c(0,1,0,1))
      X.lookup <- rmh(model=mod17,start=list(n.start=100),
                      control=list(nrep=nr,nverb=nv))
                   
   # Strauss with trend
   tr <- function(x,y){x <- x/250; y <- y/250;
   			   exp((6*x + 5*y - 18*x^2 + 12*x*y - 9*y^2)/6)
                         }
   beta <- 0.3
   gmma <- 0.5
   r    <- 45
   mod17 <- list(cif="strauss",par=c(beta=beta,gamma=gmma,r=r),w=c(0,250,0,250),
                 trend=tr3)
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
\keyword{datagen}
