\name{bw.relrisk}
\alias{bw.relrisk}
\title{
  Cross Validated Bandwidth Selection for Relative Risk Estimation
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the estimation of relative risk.
}
\usage{
bw.relrisk(X, method = "likelihood", nh = 32)
}

\arguments{
  \item{X}{
    A multitype point pattern (object of class \code{"ppp"}
    which has factor valued marks).
  }
  \item{method}{
    Character string determining the cross-validation method.
    Current options are \code{"likelihood"},
    \code{"leastsquares"} or
    \code{"weightedleastsquares"}.
  }
  \item{nh}{
    Number of trial values of smoothing bandwith \code{sigma}
    to consider.
  }
}
\details{
  This function selects an appropriate bandwidth for the nonparametric
  estimation of relative risk using \code{\link{relrisk}}.
  
  Consider the indicators \eqn{y_{ij}}{y[i,j]} which equal \eqn{1} when
  data point \eqn{x_i}{x[i]} belongs to type \eqn{j}, and equal \eqn{0}
  otherwise.
  For a particular value of smoothing bandwidth,
  let \eqn{\hat p_j(u)}{p*[j](u)} be the estimated
  probabilities that a point at location \eqn{u} will belong to
  type \eqn{j}. 
  Then the bandwidth is chosen to minimise either the likelihood,
  the squared error, or the approximately standardised squared error, of the
  indicators \eqn{y_{ij}}{y[i,j]} relative to the fitted
  values  \eqn{\hat p_j(x_i)}{p*[j](x[i])}. See Diggle (2003).

  The result is a numerical value giving the selected bandwidth \code{sigma}.
  The result also belongs to the class \code{"bw.relrisk"}
  allowing it to be printed and plotted. The plot shows the cross-validation
  criterion as a function of bandwidth.  
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.relrisk"}.
}
\seealso{
  \code{\link{relrisk}}
}
\examples{
  data(urkiola)
  b <- bw.relrisk(urkiola)
  b
  plot(b)
}
\references{
  Diggle, P.J. (2003)
  \emph{Statistical analysis of spatial point patterns},
  Second edition. Arnold.
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{methods}
\keyword{smooth}
