\name{envelope}
\alias{envelope}
\title{Simulation envelopes of summary function}
\description{
  Computes simulation envelopes of a summary function.
}
\usage{
  envelope(Y, fun=Kest, nsim=99, nrank=1, verbose=TRUE, \dots)
}
\arguments{
  \item{Y}{Either a point pattern (object of class
    \code{"ppp"}) or a fitted point process model
    (object of class \code{"ppm"}).
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
      Rank of the envelope value amongst the \code{nsim} simulated
      values. A rank of 1 means that the minimum and maximum
      simulated values will be used.
    }
    \item{verbose}{
      Logical flag indicating whether to print progress reports
      during the simulations.
    }
    \item{\dots}{
      Extra arguments passed to \code{fun}.
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
}
\details{
  Simulation envelopes can be used to assess the goodness-of-fit of
  a point process model to point pattern data. See the References.

  If \code{Y} is a point pattern (an object of class \code{"ppp"})
  then this routine generates \code{nsim} simulations of
  Complete Spatial Randomness (i.e. \code{nsim} simulated point patterns
  each being a realisation of the uniform Poisson point process)
  with the same intensity as the pattern \code{Y}.

  If \code{Y} is a fitted point process model (an object of class
  \code{"ppm"}) then this routine generates \code{nsim} simulated
  realisations of that model.
  
  The summary statistic \code{fun} is applied to each of these simulated
  patterns. Typically \code{fun} is one of the functions
  \code{Kest}, \code{Gest}, \code{Fest}, \code{Jest} or \code{pcf}.
  It can also be a home-made function; it should return an object
  of class \code{"fv"}.

  Upper and lower pointwise envelopes are computed pointwise (i.e.
  for each value of the distance argument \eqn{r}), by sorting the
  \code{nsim} simulated values, and taking the \code{m}-th lowest
  and \code{m}-th highest values, where \code{m = nrank}.
  For example if \code{nrank=1}, the upper and lower envelopes
  are the pointwise maximum and minimum of the simulated values.

  The significance level of the associated Monte Carlo test is
  \code{alpha = 2 * nrank/(1 + nsim)}. 
  
  The return value is an object of class \code{"fv"} containing
  the summary function for the data point pattern
  and the upper and lower simulation envelopes. It can be plotted
  using \code{plot.fv}.

  Arguments can be passed to the function \code{fun} through
  \code{...}. This makes it possible to select the edge correction
  used to calculate the summary statistic. See the Examples.
}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Arnold, 2003.

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
 X <- rpoispp(42)

 # Envelope of K function under CSR
 \dontrun{
 plot(envelope(X))
 }
 \testonly{
  plot(envelope(X, nsim=5))
 }

 # Translation edge correction (this is also FASTER):
 \dontrun{
 plot(envelope(X, correction="translate"))
 }
 \testonly{
  plot(envelope(X, nsim=5, correction="translate"))
 }

 # Envelope of K function for simulations from model 
 data(cells)
 fit <- ppm(cells, ~1, Strauss(0.05))
 \dontrun{
 plot(envelope(fit))
 }
 \testonly{
  plot(envelope(fit, nsim=4))
 }

 # Envelope of G function under CSR
 \dontrun{
 plot(envelope(X, Gest))
 }
 \testonly{
  plot(envelope(X, Gest, nsim=5))
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
 
 
