\name{Kest}
\alias{Kest}
\title{K-function}
\description{
Estimates the reduced second moment function \eqn{K(r)} 
from a point pattern in a window of arbitrary shape.
}
\synopsis{
  Kest(X, r=NULL, breaks=NULL, slow,
     correction=c("border", "isotropic", "Ripley", "translate"), \dots)
}
\usage{
  Kest(X)
  Kest(X, r)
  Kest(X, r, correction=c("border", "isotropic", "Ripley", "translate"))
  Kest(X, breaks=breaks)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of \eqn{K(r)} will be computed.
    An object of class \code{"ppp"}, or data
    in any format acceptable to \code{\link{as.ppp}()}.
    }
    \item{r}{
      vector of values for the argument \eqn{r} at which \eqn{K(r)} 
      should be evaluated. There is a sensible default.
    }
    \item{breaks}{
      An alternative to the argument \code{r}.
      Not normally invoked by the user.
      See Details.
    }
    \item{correction}{
      A character vector containing any selection of the
      options \code{"border"}, \code{"bord.modif"},
      \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
      It specifies the edge correction(s) to be applied.
    }
}
\value{
  A data frame containing columns
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the function \eqn{K} has been  estimated
  }
  \item{theo}{the theoretical value \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}
    for a stationary Poisson process
  }
  together with columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K(r)} obtained by the edge corrections
  named.
}
\details{
  The \eqn{K} function (variously called ``Ripley's K-function''
  and the ``reduced second moment function'')
  of a stationary point process \eqn{X} is defined so that
  \eqn{\lambda K(r)}{lambda K(r)} equals the expected number of
  additional random points within a distance \eqn{r} of a
  typical random point of \eqn{X}. Here \eqn{\lambda}{lambda}
  is the intensity of the process,
  i.e. the expected number of points of \eqn{X} per unit area.
  The \eqn{K} function is determined by the 
  second order moment properties of \eqn{X}.
 
  An estimate of \eqn{K} derived from a spatial point pattern dataset
  can be used in exploratory data analysis and formal inference
  about the pattern (Cressie, 1991; Diggle, 1983; Ripley, 1988).
  In exploratory analyses, the estimate of \eqn{K} is a useful statistic 
  summarising aspects of inter-point ``dependence'' and ``clustering''.
  For inferential purposes, the estimate of \eqn{K} is usually compared to the 
  true value of \eqn{K} for a completely random (Poisson) point process,
  which is \eqn{K(r) = \pi r^2}{K(r) = pi * r^2}.
  Deviations between the empirical and theoretical \eqn{K} curves
  may suggest spatial clustering or spatial regularity.
 
  This routine \code{Kest} estimates the \eqn{K} function
  of a stationary point process, given observation of the process
  inside a known, bounded window. 
  The argument \code{X} is interpreted as a point pattern object 
  (of class \code{"ppp"}, see \code{\link{ppp.object}}) and can
  be supplied in any of the formats recognised by
  \code{\link{as.ppp}()}.

  The estimation of \eqn{K} is hampered by edge effects arising from 
  the unobservability of points of the random pattern outside the window. 
  An edge correction is needed to reduce bias (Baddeley, 1998; Ripley, 1988). 
  The corrections implemented here are
  \describe{
    \item{border}{the border method or
      ``reduced sample'' estimator (see Ripley, 1988). This is
      the least efficient (statistically) and the fastest to compute.
      It can be computed for a window of arbitrary shape.
    }
    \item{isotropic/Ripley}{Ripley's isotropic correction
      (see Ripley, 1988; Ohser, 1983).
      This is currently implemented only for rectangular windows.
    }
    \item{translate}{Translation correction (Ohser, 1983).
      Implemented for all window geometries, but slow for
      complex windows. 
    }
  }
  Note that the estimator assumes the process is stationary (spatially
  homogeneous). For inhomogeneous point patterns, see
  \code{\link{Kinhom}}.

  The estimator \code{Kest} ignores marks.
  Its counterparts for multitype point patterns
  are \code{\link{Kcross}}, \code{\link{Kdot}},
  and for general marked point patterns
  see \code{\link{Kmulti}}. 

  Some writers, particularly Stoyan (1994, 1995) advocate the use of
  the ``pair correlation function''
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}.
  See \code{\link{pcf}} on how to estimate this function.

}
\references
Baddeley, A.J. Spatial sampling and censoring.
     In O.E. Barndorff-Nielsen, W.S. Kendall and
     M.N.M. van Lieshout (eds) 
     \emph{Stochastic Geometry: Likelihood and Computation}.
     Chapman and Hall, 1998.
     Chapter 2, pages 37--78.
  
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
  Academic Press, 1983.

  Ohser, J. (1983)
  On estimators for the reduced second moment measure of
  point processes. \emph{Mathematische Operationsforschung und
  Statistik, series Statistics}, \bold{14}, 63 -- 71.
    
  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
} 
\section{Warnings}{
  The estimator of \eqn{K(r)} is approximately unbiased for each fixed \eqn{r}.
  Bias increases with \eqn{r} and depends on the window geometry.
  For a rectangular window it is prudent to restrict the \eqn{r} values to
  a maximum of \eqn{1/4} of the smaller side length of the rectangle.
  Bias may become appreciable for point patterns consisting of 
  fewer than 15 points.
 
  While \eqn{K(r)} is always a non-decreasing function, the estimator 
  of \eqn{K} is not guaranteed to be non-decreasing. This is rarely 
  a problem in practice.
}
\seealso{
  \code{\link{Fest}},
  \code{\link{Gest}},
  \code{\link{Jest}},
  \code{\link{pcf}},
  \code{\link{reduced.sample}},
  \code{\link{Kcross}},
  \code{\link{Kdot}},
  \code{\link{Kinhom}},
  \code{\link{Kmulti}}
}
\examples{
 library(spatstat)
 pp <- runifpoint(50)
 Kpp <- Kest(pp)
 \dontrun{
 plot(Kpp$r, Kpp$border, type="l", xlab="r", ylab="K(r)", ylim=c(0,1),
       main = "K-function")
 r <- Kpp$r
 lines(r, Kpp$theo, lty=2)
 legend(0.5, 2, c("empirical (border correction)", "Poisson process"), lty=c(1,2))
 }
 data(cells)
 K <- Kest(cells, correction="isotropic")
 \dontrun{
 plot(K$r, K$iso, type="l")
 plot(iso ~ r, type="l", data=K)
 # restrict the plot to values of r less than 0.1
 plot(iso ~ r, type="l", data=K[K$r <= 0.1, ])
 plot(iso ~ r, type="l", data=K, subset=(r <= 0.1))
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
 
 
