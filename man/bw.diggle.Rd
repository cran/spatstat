\name{bw.diggle}
\alias{bw.diggle}
\title{
  Cross Validated Bandwidth Selection for Kernel Density
}
\description{
  Uses cross-validation to select a smoothing bandwidth
  for the kernel estimation of point process intensity.
}
\usage{
   bw.diggle(X)
}
\arguments{
  \item{X}{
    A point pattern (object of class \code{"ppp"}).
  }
}
\details{
  This function selects an appropriate bandwidth \code{sigma}
  for the kernel estimator of point process intensity
  computed by \code{\link{density.ppp}}.

  The bandwidth \eqn{\sigma}{sigma} is chosen to 
  minimise the mean-square error criterion defined by Diggle (1985).
  The algorithm computes the mean-square error by the method
  of Berman and Diggle (1989). See Diggle (2003, page 117) for a
  summary.

  The result is a numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted to show the (rescaled) mean-square error
  as a function of \code{sigma}.
}
\value{
  A numerical value giving the selected bandwidth.
  The result also belongs to the class \code{"bw.optim"}
  which can be plotted.
}
\seealso{
  \code{\link{density.ppp}}
}
\examples{
  data(lansing)
  attach(split(lansing))
  b <- bw.diggle(hickory)
  plot(b, ylim=c(-2.5, 0), main="Cross validation for hickories")
  plot(density(hickory, b))
}
\references{
  Berman, M. and Diggle, P. (1989)
  Estimating weighted integrals of the
  second-order intensity of a spatial point process.
  \emph{Journal of the Royal Statistical Society, series B}
  \bold{51}, 81--92.

  Diggle, P.J. (1985)
  A kernel method for smoothing point process data.
  \emph{Applied Statistics} (Journal of the Royal Statistical Society,
  Series C) \bold{34} (1985) 138--147.

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
