\name{F3est}
\Rdversion{1.1}
\alias{F3est}
\title{
  Empty Space Function of a Three-Dimensional Point Pattern
}
\description{
  Estimates the empty space function \eqn{F_3(r)}{F3(r)} from 
  a three-dimensional point pattern.
}
\usage{
F3est(X, ..., rmax = NULL, nrval = 128, vside = NULL, correction = c("rs", "km", "cs"))
}
\arguments{
  \item{X}{
    Three-dimensional point pattern (object of class \code{"pp3"}).
  }
  \item{\dots}{
    Ignored.
  }
  \item{rmax}{
    Optional. Maximum value of argument \eqn{r} for which
    \eqn{F_3(r)}{F3(r)} will be estimated. 
  }
  \item{nrval}{
    Optional. Number of values of \eqn{r} for which
    \eqn{F_3(r)}{F3(r)} will be estimated. A large value of \code{nrval}
    is required to avoid discretisation effects.
  }
  \item{vside}{
    Optional. 
    Side length of the voxels in the discrete approximation.
  }
  \item{correction}{
    Optional. Character vector specifying the edge correction(s)
    to be applied. See Details.
  }
}
\details{
  For a stationary point process \eqn{\Phi}{Phi} in three-dimensional
  space, the empty space function is
  \deqn{
    F_3(r) = P(d(0,\Phi) \le r)
  }{
    F3(r) = P(d(0,Phi) <= r)
  }
  where \eqn{d(0,\Phi)}{d(0,Phi)} denotes the distance from a fixed
  origin \eqn{0} to the nearest point of \eqn{\Phi}{Phi}.
  
  The three-dimensional point pattern \code{X} is assumed to be a
  partial realisation of a stationary point process \eqn{\Phi}{Phi}.
  The empty space function of \eqn{\Phi}{Phi} can then be estimated using
  techniques described in the References.

  The box containing the point
  pattern is discretised into cubic voxels of side length \code{vside}.
  The distance function \eqn{d(u,\Phi)}{d(u,Phi)} is computed for
  every voxel centre point
  \eqn{u} using a three-dimensional version of the distance transform
  algorithm (Borgefors, 1986). The empirical cumulative distribution
  function of these values, with appropriate edge corrections, is the
  estimate of \eqn{F_3(r)}{F3(r)}.

  The available edge corrections are:
  \describe{
    \item{\code{"rs"}:}{
      the reduced sample (aka minus sampling, border correction)
      estimator (Baddeley et al, 1993)
    }
    \item{\code{"km"}:}{
      the three-dimensional version of the
      Kaplan-Meier estimator (Baddeley and Gill, 1997)
    }
    \item{\code{"cs"}:}{
      the three-dimensional generalisation of
      the Chiu-Stoyan or Hanisch estimator (Chiu and Stoyan, 1998).
    }
  }
}
\value{
  A function value table (object of class \code{"fv"}) that can be
  plotted, printed or coerced to a data frame containing the function values.
}
\references{
  Baddeley, A.J, Moyeed, R.A., Howard, C.V. and Boyde, A.
  Analysis of a three-dimensional point pattern with replication.
  \emph{Applied Statistics} \bold{42} (1993) 641--668.

  Baddeley, A.J. and Gill, R.D. (1997)
  Kaplan-Meier estimators of interpoint distance
  distributions for spatial point processes.
  \emph{Annals of Statistics} \bold{25}, 263--292.

  Borgefors, G. (1986)
  Distance transformations in digital images.
  \emph{Computer Vision, Graphics and Image Processing}
  \bold{34}, 344--371.

  Chiu, S.N. and Stoyan, D. (1998)
  Estimators of distance distributions for spatial patterns.
  \emph{Statistica Neerlandica} \bold{52}, 239--246.
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rana Moyeed.
}
\section{Warnings}{
  A large value of \code{nrval} is required in order to avoid
  discretisation effects (due to the use of histograms in the
  calculation).
}
\seealso{
  \code{\link{G3est}},
  \code{\link{K3est}}.
}
\examples{
  X <- pp3(runif(42), runif(42), runif(42), box3(c(0,1)))
  plot(F3est(X))
}
\keyword{spatial}
\keyword{nonparametric}
