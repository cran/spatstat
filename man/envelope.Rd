\name{envelope}
\alias{envelope}
\title{Simulation envelopes of summary function}
\description{
  Computes simulation envelopes of a summary function.
}
\usage{
  envelope(Y, fun=Kest, nsim=99, nrank=1, \dots, 
  simulate=NULL, verbose=TRUE, clipdata=TRUE,
  start=NULL,control=list(nrep=1e5,expand=1.5),
  transform=NULL,global=FALSE,ginterval=NULL,
  savefuns=FALSE, savepatterns=FALSE,
  nsim2=nsim, VARIANCE=FALSE, nSD=2, Yname=NULL, internal=NULL)
}
\arguments{
  \item{Y}{Either a point pattern (object of class
    \code{"ppp"}) or a fitted point process model
    (object of class \code{"ppm"} or \code{"kppm"}).
    }
    \item{fun}{
      Function that computes the desired summary statistic
      for a point pattern. 
    }
    \item{nsim}{
      Number of simulated point patterns to be generated
      when computing the envelopes.
    }
    \item{nrank}{
      Integer. Rank of the envelope value amongst the \code{nsim} simulated
      values. A rank of 1 means that the minimum and maximum
      simulated values will be used.
    }
    \item{\dots}{
      Extra arguments passed to \code{fun}.
    }
    \item{simulate}{
      Optional. Specifies how to generate the simulated point patterns.
      If \code{simulate} is an expression in the R language, then this
      expression will be evaluated \code{nsim} times,
      to obtain \code{nsim} point patterns which are taken as the
      simulated patterns from which the envelopes are computed.
      If \code{simulate} is a list of point patterns, then the entries
      in this list will be treated as the simulated patterns from which
      the envelopes are computed.
      Alternatively \code{simulate} may be an object produced by the
      \code{envelope} command: see Details.
    }
    \item{verbose}{
      Logical flag indicating whether to print progress reports
      during the simulations.
    }
    \item{clipdata}{
      Logical flag indicating whether the data point pattern should be
      clipped to the same window as the simulated patterns,
      before the summary function for the data is computed.
      This should usually be \code{TRUE} to ensure that the
      data and simulations are properly comparable.
    }
    \item{start,control}{
      Optional. These specify the arguments \code{start} and \code{control}
      of \code{rmh}, giving complete control over the
      simulation algorithm.
    }
    \item{transform}{
      Optional. A transformation to be applied to the
      function values, before the envelopes are computed.
      An expression object (see Details).
    }
    \item{global}{
      Logical flag indicating whether envelopes should be pointwise
      (\code{global=FALSE}) or simultaneous (\code{global=TRUE}).
    }
    \item{ginterval}{
      Optional.
      A vector of length 2 specifying
      the interval of \eqn{r} values for the simultaneous critical
      envelopes. Only relevant if \code{global=TRUE}.
    }
    \item{savefuns}{
      Logical flag indicating whether to save all the simulated
      function values.
    }
    \item{savepatterns}{
      Logical flag indicating whether to save all the simulated
      point patterns.
    }
    \item{nsim2}{
      Number of extra simulated point patterns to be generated
      if it is necessary to use simulation to estimate the theoretical
      mean of the summary function. Only relevant when \code{global=TRUE}
      and the simulations are not based on CSR.
    }
    \item{VARIANCE}{
      Logical. If \code{TRUE}, critical envelopes will be calculated
      as sample mean plus or minus \code{nSD} times sample standard
      deviation.
    }
    \item{nSD}{
      Number of estimated standard deviations used to determine
      the critical envelopes, if \code{VARIANCE=TRUE}.
    }
    \item{Yname}{
      Character string that should be used as the name of the 
      data point pattern \code{Y} when printing or plotting the results.
    }
    \item{internal}{
      Do not use this argument.
      It is used by the package's internal code.
    }
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the summary function \code{fun} has been  estimated
  }
  \item{obs}{
    values of the summary function for the data point pattern
  }
  \item{lo}{
    lower envelope of simulations
  }
  \item{hi}{
    upper envelope of simulations
  }
  and \emph{either}
  \item{theo}{
    theoretical value of the summary function under CSR
    (Complete Spatial Randomness, a uniform Poisson point process)
    if the simulations were generated according to CSR
  }
  \item{mmean}{
    estimated theoretical value of the summary function,
    computed by averaging simulated values, 
    if the simulations were not generated according to CSR.
  }
  Additionally, if \code{savepatterns=TRUE}, the return value has an attribute
  \code{"simpatterns"} which is a list containing the \code{nsim}
  simulated patterns. If \code{savefuns=TRUE}, the return value
  has an attribute \code{"simfuns"} which is an object of class
  \code{"fv"} containing the summary functions
  computed for each of the \code{nsim} simulated patterns.
}
\details{
  Simulation envelopes can be used to assess the goodness-of-fit of
  a point process model to point pattern data. See the References.

  This function first generates \code{nsim} random point patterns
  in one of the following ways. 
  \itemize{
    \item 
    If \code{Y} is a point pattern (an object of class \code{"ppp"})
    and \code{simulate=NULL},
    then this routine generates \code{nsim} simulations of
    Complete Spatial Randomness (i.e. \code{nsim} simulated point patterns
    each being a realisation of the uniform Poisson point process)
    with the same intensity as the pattern \code{Y}.
    (If \code{Y} is a multitype point pattern, then the simulated patterns
    are also given independent random marks; the probability
    distribution of the random marks is determined by the
    relative frequencies of marks in \code{Y}.) 
    \item
    If \code{Y} is a fitted point process model (an object of class
    \code{"ppm"} or \code{"kppm"}) and \code{simulate=NULL},
    then this routine generates \code{nsim} simulated
    realisations of that model.
    \item
    If \code{simulate} is supplied, then it determines how the
    simulated point patterns are generated. It may be either
    \itemize{
      \item
      an expression in the R language, typically containing a call
      to a random generator. This expression will be evaluated
      \code{nsim} times to yield \code{nsim} point patterns. For example
      if \code{simulate=expression(runifpoint(100))} then each simulated
      pattern consists of exactly 100 independent uniform random points.
      \item
      a list of point patterns.
      The entries in this list will be taken as the simulated patterns.
      \item
      an object of class \code{"envelope"}. This should have been
      produced by calling \code{envelope} with the
      argument \code{savepatterns=TRUE}.
      The simulated point patterns that were saved in this object
      will be extracted and used as the simulated patterns for the
      new envelope computation. This makes it possible to plot envelopes
      for two different summary functions based on exactly the same set of
      simulated point patterns.
      }
  }
  
  The summary statistic \code{fun} is applied to each of these simulated
  patterns. Typically \code{fun} is one of the functions
  \code{Kest}, \code{Gest}, \code{Fest}, \code{Jest}, \code{pcf},
  \code{Kcross}, \code{Kdot}, \code{Gcross}, \code{Gdot},
  \code{Jcross}, \code{Jdot}, \code{Kmulti}, \code{Gmulti},
  \code{Jmulti} or \code{Kinhom}. It may also be a character string
  containing the name of one of these functions.

  The statistic \code{fun} can also be a user-supplied function;
  if so, then it must have arguments \code{X} and \code{r}
  like those in the functions listed above, and it must return an object
  of class \code{"fv"}.

  Upper and lower critical envelopes are computed in one of the following ways:
  \describe{
    \item{pointwise:}{by default, envelopes are calculated pointwise
      (i.e. for each value of the distance argument \eqn{r}), by sorting the
      \code{nsim} simulated values, and taking the \code{m}-th lowest
      and \code{m}-th highest values, where \code{m = nrank}.
      For example if \code{nrank=1}, the upper and lower envelopes
      are the pointwise maximum and minimum of the simulated values.

      The pointwise envelopes are \bold{not} \dQuote{confidence bands}
      for the true value of the function! Rather,
      they specify the critical points for a Monte Carlo test
      (Ripley, 1981). The test is constructed by choosing a
      \emph{fixed} value of \eqn{r}, and rejecting the null hypothesis if the
      observed function value
      lies outside the envelope \emph{at this value of} \eqn{r}.
      This test has exact significance level
      \code{alpha = 2 * nrank/(1 + nsim)}.
    }
    \item{simultaneous:}{if \code{global=TRUE}, then the envelopes are
      determined as follows. First we calculate the theoretical mean value of
      the summary statistic (if we are testing CSR, the theoretical
      value is supplied by \code{fun}; otherwise we perform a separate
      set of \code{nsim2} simulations, compute the
      average of all these simulated values, and take this average
      as an estimate of the theoretical mean value). Then, for each simulation,
      we compare the simulated curve to the theoretical curve, and compute the
      maximum absolute difference between them (over the interval
      of \eqn{r} values specified by \code{ginterval}). This gives a
      deviation value \eqn{d_i}{d[i]} for each of the \code{nsim}
      simulations. Finally we take the \code{m}-th largest of the
      deviation values, where \code{m=nrank}, and call this
      \code{dcrit}. Then the simultaneous envelopes are of the form
      \code{lo = expected - dcrit} and \code{hi = expected + dcrit} where
      \code{expected} is either the theoretical mean value \code{theo}
      (if we are testing CSR) or the estimated theoretical value
      \code{mmean} (if we are testing another model). The simultaneous critical
      envelopes have constant width \code{2 * dcrit}.

      The simultaneous critical envelopes allow us to perform a different
      Monte Carlo test (Ripley, 1981). The test rejects the null
      hypothesis if the graph of the observed function
      lies outside the envelope \bold{at any value of} \eqn{r}.
      This test has exact significance level
      \code{alpha = nrank/(1 + nsim)}.
    }
    \item{based on sample moments:}{if \code{VARIANCE=TRUE},
      the algorithm calculates the
      (pointwise) sample mean and sample variance of
      the simulated functions. Then the envelopes are computed
      as mean plus or minus \code{nSD} standard deviations.
      These envelopes do not have an exact significance interpretation.
      They are a naive approximation to
      the critical points of the Neyman-Pearson test
      assuming the summary statistic is approximately Normally
      distributed.
    }
  }
  
  The return value is an object of class \code{"fv"} containing
  the summary function for the data point pattern,
  the upper and lower simulation envelopes, and 
  the theoretical expected value (exact or estimated) of the summary function 
  for the model being tested. It can be plotted using \code{plot.fv}.

  If \code{VARIANCE=TRUE} then the return value also includes the
  sample mean, sample variance and other quantities.

  Arguments can be passed to the function \code{fun} through
  \code{...}. This makes it possible to select the edge correction
  used to calculate the summary statistic. See the Examples.
  Selecting only a single edge
  correction will make the code run much faster.

  If \code{Y} is a fitted cluster point process model (object of
  class \code{"kppm"}), and \code{simulate=NULL},
  then the model is simulated directly
  using \code{\link{simulate.kppm}}.
  
  If \code{Y} is a fitted Gibbs point process model (object of
  class \code{"ppm"}), and \code{simulate=NULL},
  then the model is simulated
  by running the Metropolis-Hastings algorithm \code{\link{rmh}}.
  Complete control over this algorithm is provided by the 
  arguments \code{start} and \code{control} which are passed
  to \code{\link{rmh}}.

  For simultaneous critical envelopes (\code{global=TRUE})
  the following options are also useful:
  \describe{
    \item{\code{ginterval}}{determines the interval of \eqn{r} values
      over which the deviation between curves is calculated.
      It should be a numeric vector of length 2.
      There is a sensible default (namely, the recommended plotting
      interval for \code{fun(X)}, or the range of \code{r} values if
      \code{r} is explicitly specified).
    }
    \item{\code{transform}}{specifies a transformation of the
      summary function \code{fun} that will be carried out before the
      deviations are computed. It must be an expression object
      using the symbol \code{.} to represent the function value.
      For example, 
      the conventional way to normalise the \eqn{K} function
      (Ripley, 1981) is to transform it to the \eqn{L} function
      \eqn{L(r) = \sqrt{K(r)/\pi}}{L(r) = sqrt(K(r)/pi)}
      and this is implemented by setting
      \code{transform=expression(sqrt(./pi))}.
      Such transforms are only useful if \code{global=TRUE}.
    }
  }

  It is also possible to extract the summary functions for each of the
  individual simulated point patterns, by setting \code{savefuns=TRUE}.
  Then the return value also 
  has an attribute \code{"simfuns"} containing all the 
  summary functions for the individual simulated patterns.
  It is an \code{"fv"} object containing
  functions named \code{sim1, sim2, ...} representing the \code{nsim}
  summary functions.

  It is also possible to save the simulated point patterns themselves,
  by setting \code{savepatterns=TRUE}. Then the return value also has
  an attribute \code{"simpatterns"} which is a list of length
  \code{nsim} containing all the simulated point patterns.
}
\section{Warning}{
  An error may be generated if one of the simulations produces a
  point pattern that is empty, or is otherwise unacceptable to the
  function \code{fun}.
}
  
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Arnold, 2003.

  Ripley, B.D. (1981) 
  \emph{Spatial statistics}.
  John Wiley and Sons.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\seealso{
  \code{\link{fv.object}},
  \code{\link{plot.fv}},
  \code{\link{Kest}},
  \code{\link{Gest}},
  \code{\link{Fest}},
  \code{\link{Jest}},
  \code{\link{pcf}},
  \code{\link{ppp}},
  \code{\link{ppm}}
}
\examples{
 data(simdat)
 X <- simdat

 # Envelope of K function under CSR
 \dontrun{
 plot(envelope(X))
 }
 \testonly{
  plot(envelope(X, nsim=4))
 }

 # Translation edge correction (this is also FASTER):
 \dontrun{
 plot(envelope(X, correction="translate"))
 }
 \testonly{
  plot(envelope(X, nsim=4, correction="translate"))
 }

 # Envelope of K function for simulations from Gibbs model 
 data(cells)
 fit <- ppm(cells, ~1, Strauss(0.05))
 \dontrun{
 plot(envelope(fit))
 plot(envelope(fit), global=TRUE)
 }
 \testonly{
  plot(envelope(fit, nsim=4))
  plot(envelope(fit, nsim=4, global=TRUE))
 }

 # Envelope of K function for simulations from cluster model 
 data(redwood)
 fit <- kppm(redwood, ~1, "Thomas")
 \dontrun{
 plot(envelope(fit, Gest))
 plot(envelope(fit, Gest, global=TRUE))
 }
 \testonly{
  plot(envelope(fit, Gest, nsim=4))
  plot(envelope(fit, Gest, nsim=4, global=TRUE))
 }

 # Envelope of G function under CSR
 \dontrun{
 plot(envelope(X, Gest))
 }
 \testonly{
  plot(envelope(X, Gest, nsim=4))
 }

 # Envelope of L function under CSR
 #  L(r) = sqrt(K(r)/pi)
 \dontrun{
  E <- envelope(X, Kest)
  plot(E, sqrt(./pi) ~ r)
 }
 \testonly{
  E <- envelope(X, Kest, nsim=4)
  plot(E, sqrt(./pi) ~ r)
 }

 # Simultaneous critical envelope for L function
 # (alternatively, use Lest)
 \dontrun{
  plot(envelope(X, Kest, transform=expression(sqrt(./pi)), global=TRUE))
 }
 \testonly{
 plot(envelope(X, Kest, nsim=4,transform=expression(sqrt(./pi)), global=TRUE))
 }

 # How to pass arguments needed to compute the summary functions:
 # We want envelopes for Jcross(X, "A", "B") 
 # where "A" and "B" are types of points in the dataset 'demopat'

 data(demopat)
 \dontrun{
 plot(envelope(demopat, Jcross, i="A", j="B"))
 }
 \testonly{
 plot(envelope(demopat, Jcross, i="A", j="B", nsim=4))
 }
 
 # Use of `simulate'
 \dontrun{
 plot(envelope(cells, Gest, simulate=expression(runifpoint(42))))
 plot(envelope(cells, Gest, simulate=expression(rMaternI(100,0.02))))
 }
 \testonly{
  plot(envelope(cells, Gest, simulate=expression(runifpoint(42)), nsim=4))
  plot(envelope(cells, Gest, simulate=expression(rMaternI(100, 0.02)), nsim=4))
  plot(envelope(cells, Gest, simulate=expression(runifpoint(42)),
nsim=4, global=TRUE))
  plot(envelope(cells, Gest, simulate=expression(rMaternI(100, 0.02)),
nsim=4, global=TRUE))
 }

 # Envelope under random toroidal shifts
 data(amacrine)
 \dontrun{
 plot(envelope(amacrine, Kcross, i="on", j="off",
               simulate=expression(rshift(amacrine, radius=0.25)))) 
 }

 # Envelope under random shifts with erosion
 \dontrun{
 plot(envelope(amacrine, Kcross, i="on", j="off",
              simulate=expression(rshift(amacrine, radius=0.1, edge="erode"))))
 }
  
 # Envelope of INHOMOGENEOUS K-function with fitted trend
\dontrun{
 trend <- density.ppp(X, 1.5)
 plot(envelope(X, Kinhom, lambda=trend,
         simulate=expression(rpoispp(trend))))
 }

 # Precomputed list of point patterns
 X <- rpoispp(50)
 PatList <- list()
 for(i in 1:20) PatList[[i]] <- runifpoint(X$n)
 plot(envelope(X, Kest, nsim=20, simulate=PatList))

# re-using the same point patterns
 EK <- envelope(X, Kest, nsim=10, savepatterns=TRUE)
 EG <- envelope(X, Kest, nsim=10, simulate=EK)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{htest}
\keyword{hplot}
\keyword{iteration}
 
 
