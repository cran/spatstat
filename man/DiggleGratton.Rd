\name{DiggleGratton}
\alias{DiggleGratton}
\title{Diggle-Gratton model}
\description{
  Creates an instance of the Diggle-Gratton pairwise interaction
  point process model, which can then be fitted to point pattern data.
}
\usage{
  DiggleGratton(delta, rho)
}
\arguments{
  \item{delta}{lower threshold \eqn{\delta}{delta}}
  \item{rho}{upper threshold \eqn{\rho}{rho}}
}
\value{
  An object of class \code{"interact"}
  describing the interpoint interaction
  structure of a point process. 
}
\details{
  Diggle and Gratton (1984, pages 208-210)
  introduced the pairwise interaction point
  process with pair potential \eqn{h(t)} of the form
  \deqn{
    h(t) = \left( \frac{t-\delta}{\rho-\delta} \right)^\beta
    \quad\quad \mbox{  if  } \delta \le t \le \rho
  }{
    h(t) = ((t - delta)/(rho - delta))^beta, {    } delta <= t <= rho
  }
  with \eqn{h(t) = 0} for \eqn{t < \delta}{t < delta}
  and  \eqn{h(t) = 1} for \eqn{t > \rho}{t > rho}.

  Here \eqn{\delta \le \rho}{delta <= rho} are irregular parameters.
  The potential is inhibitory, i.e.\ this model is only appropriate for
  regular point patterns.

  Note that the irregular parameters 
  \eqn{\delta, \rho}{delta, rho} must be fixed, while the
  regular parameter \eqn{\beta}{beta} will be estimated. 
}
\seealso{
  \code{\link{mpl}},
  \code{\link{ppm.object}},
  \code{\link{Pairwise}}
}
\examples{
   data(cells)
   mpl(cells, ~1, DiggleGratton(0.05, 0.1))
   # NOTE: resulting estimate of beta is 8.569, not exp(8.569) = 5267.54
}
\references{
  Diggle, P.J. and Gratton, R.J. (1984)
  Monte Carlo methods of inference for implicit statistical models.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{46}, 193 -- 212.
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}

