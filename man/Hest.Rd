\name{Hest}
\alias{Hest}
\title{Spherical Contact Distribution Function}
\description{
  Estimates the spherical contact distribution function of a
  random set.
}
\usage{
Hest(X, ...)
}
\arguments{
  \item{X}{The observed random set.
    An object of class \code{"ppp"}, \code{"psp"} or \code{"owin"}.
  }
  \item{\dots}{Arguments passed to \code{\link{as.mask}}
    to control the discretisation.}
}
\details{
  The spherical contact distribution function
  of a stationary random set \eqn{X}
  is the cumulative distribution function \eqn{H} of the distance
  from a fixed point in space to the nearest point of \eqn{X},
  given that the point lies outside \eqn{X}.
  That is, \eqn{H(r)} equals
  the probability that \code{X} lies closer than \eqn{r} units away
  from the fixed point \eqn{x}, given that \code{X} does not cover \eqn{x}.

  For a point process, the spherical contact distribution function
  is the same as the empty space function \eqn{F} discussed
  in \code{\link{Fest}}. 

  For \code{Hest}, the argument \code{X} may be a point pattern
  (object of class \code{"ppp"}), a line segment pattern
  (object of class \code{"psp"}) or a window (object of class
  \code{"owin"}). It is assumed to be a realisation of a stationary
  random set.

  The algorithm first calls \code{\link{distmap}} to compute the
  distance transform of \code{X}, then computes the Kaplan-Meier
  and reduced-sample estimates of the cumulative distribution
  following Hansen et al (1999).
}
\value{
  An object of class \code{"fv"}, see \code{\link{fv.object}},
  which can be plotted directly using \code{\link{plot.fv}}.

  Essentially a data frame containing five columns:
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{H(r)} has been  estimated
  }
  \item{rs}{the ``reduced sample'' or ``border correction''
    estimator of \eqn{H(r)}
  }
  \item{km}{the spatial Kaplan-Meier estimator of \eqn{H(r)}
  }
  \item{hazard}{the hazard rate \eqn{\lambda(r)}{lambda(r)}
    of \eqn{H(r)} by the spatial Kaplan-Meier method
  }
  \item{raw}{the uncorrected estimate of \eqn{H(r)},
  i.e. the empirical distribution of the distance from 
  a fixed point in the window to the nearest point of \code{X}
  }
}
\references{
  Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37-78.
  
  Baddeley, A.J. and Gill, R.D. 
    The empty space hazard of a spatial pattern.
    Research Report 1994/3, Department of Mathematics,
    University of Western Australia, May 1994.

  Hansen, M.B., Baddeley, A.J. and Gill, R.D.
  First contact distributions for spatial patterns:
  regularity and estimation.
  \emph{Advances in Applied Probability} \bold{31} (1999) 15-33.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.
}
\seealso{\code{\link{Fest}}}
\examples{
   X <- runifpoint(42)
   H <- Hest(X)
   Y <- rpoisline(10)
   H <- Hest(Y)
   data(heather)
   H <- Hest(heather$medium)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
