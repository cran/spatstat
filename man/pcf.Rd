\name{pcf}
\alias{pcf}
\title{Pair Correlation Function}
\description{
  Estimates the pair correlation function of
  a point pattern.
}
\usage{
 pcf(X, \dots, method="c")
}
\arguments{
  \item{X}{
    Either the observed data point pattern,
    or an estimate of its \eqn{K} function,
    or an array of multitype \eqn{K} functions
    (see Details).
  }
  \item{\dots}{
    Arguments controlling the smoothing spline
    function \code{smooth.spline}.
  }
  \item{method}{
    Letter \code{"a"}, \code{"b"}, code{"c"} or \code{"d"} indicating the
    method for deriving the pair correlation function from the
    \code{K} function.
  }
} 
\value{
  Either a function value table
  (object of class \code{"fv"}, see \code{\link{fv.object}})
  representing a pair correlation function,
  or a function array (object of class \code{"fasp"},
  see \code{\link{fasp.object}})
  representing an array of pair correlation functions.

  If \code{X} is an unmarked point pattern, the return value is
  a function value table (class \code{"fv"}). It is
  essentially a data frame containing (at least) the variables
  \item{r}{the vector of values of the argument \eqn{r} 
    at which the pair correlation function \eqn{g(r)} has been  estimated
  }
  \item{pcf}{vector of values of \eqn{g(r)}
  }
  If code{X} is a function value table (class \code{"fv"}) representing
  the estimated \eqn{K} function of a point pattern
  (obtained from \code{\link{Kest}}, \code{\link{Kinhom}}
  or \code{\link{Kest.fft}}),
  then the return value is again a function value table,
  representing the pair correlation function.

  If \code{X} is a multitype point pattern,
  the return value is a function array (class \code{"fasp"}).
  This can be thought of as a matrix \code{Y} each of whose entries
  \code{Y[i,j]} is a function value table (class \code{"fv"})
  representing the pair correlation function between
  points of type \code{i} and points of type \code{j}.

  If \code{X} is a function array (class \code{"fasp"})
  representing an array of estimated \eqn{K} functions
  (obtained from \code{\link{alltypes}}),
  the return value is another function array, containing
  the corresponding pair correlation functions.
}
\details{
  The pair correlation function of a stationary point process is
  \deqn{
    g(r) = \frac{K'(r)}{2\pi r}
  }{
    g(r) = K'(r)/ ( 2 * pi * r) 
  }
  where \eqn{K'(r)} is the derivative of \eqn{K(r)}, the
  reduced second moment function (aka ``Ripley's \eqn{K} function'')
  of the point process. See \code{\link{Kest}} for information
  about \eqn{K(r)}. For a stationary Poisson process, the
  pair correlation function is identically equal to 1. Values
  \eqn{g(r) < 1} suggest inhibition between points;
  values greater than 1 suggest clustering.

  We also apply the same definition to
  other variants of the classical \eqn{K} function,
  such as the multitype \eqn{K} functions
  (see \code{\link{Kcross}}, \code{\link{Kdot}}) and the
  inhomogeneous \eqn{K} function (see \code{\link{Kinhom}}).
  For all these variants, the benchmark value of
  \eqn{K(r) = \pi r^2}{K(r) = pi * r^2} corresponds to
  \eqn{g(r) = 1}.

  This routine computes an estimate of \eqn{g(r)}
  from an estimate of \eqn{K(r)} or its variants,
  using smoothing splines to approximate the derivative.

  The argument \code{X} may be either
  \itemize{
    \item
    a point pattern for which an estimate of the pair correlation function
    should be computed. This should be an object of class \code{"ppp"},
    or in a format recognised by \code{\link{as.ppp}()}.
    \item
    an estimated \eqn{K} function,
    given as a function value table (object of class \code{"fv"},
    see \code{\link{fv.object}}).
    This object should be the value returned by
    \code{\link{Kest}}, \code{\link{Kcross}}, \code{\link{Kmulti}}
    or \code{\link{Kinhom}}.
    \item
    a function array (object of class \code{"fasp"},
    see \code{\link{fasp.object}})
    containing several estimates of \eqn{K} functions.
    This should have been obtained from \code{\link{alltypes}}
    with the argument \code{fun="K"}.
  }
  If \code{X} is a point pattern, the \eqn{K} function is
  first estimated by \code{\link{Kest}}.
  
  The smoothing spline operations are performed by
  \code{\link{smooth.spline}} and \code{\link{predict.smooth.spline}}
  from the \code{modreg} library.
  Four numerical methods are available:
  \itemize{
    \item
    \bold{"a"} apply smoothing to \eqn{K(r)},
    estimate its derivative, and plug in to the formula above;
    \item 
    \bold{"b"} apply smoothing to
    \eqn{Y(r) = \frac{K(r)}{2 \pi r}}{Y(r) = K(r)/(2 * pi * r)}
    constraining \eqn{Y(0) = 0},
    estimate the derivative of \eqn{Y}, and solve;
    \item
    \bold{"c"} apply smoothing to 
    \eqn{Z(r) = \frac{K(r)}{\pi r^2}}{Y(r) = K(r)/(pi * r^2)}
    constraining \eqn{Z(0)=1},
    estimate its derivative, and solve.
    \item
    \bold{"d"} apply smoothing to 
    \eqn{V(r) = \sqrt{K(r)}}{V(r) = sqrt(K(r))},
    estimate its derivative, and solve.
  }
  Method \code{"c"} seems to be the best at 
  suppressing variability for small values of \eqn{r}.
  However it effectively constrains \eqn{g(0) = 1}.
  If the point pattern seems to have inhibition at small distances,
  you may wish to experiment with method \code{"b"} which effectively
  constrains \eqn{g(0)=0}. Method \code{"a"} seems
  comparatively unreliable.

  Useful arguments to control the splines
  include the smoothing tradeoff parameter \code{spar}
  and the degrees of freedom \code{df}. See \code{\link{smooth.spline}}
  for details.
}
\references{
  Stoyan, D, Kendall, W.S. and Mecke, J. (1995)
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag.

  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{Kcross}},
  \code{\link{Kdot}},
  \code{\link{Kmulti}},
  \code{\link{alltypes}},
  \code{\link{smooth.spline}},
  \code{\link{predict.smooth.spline}}
}
\examples{
  # univariate point pattern
  data(simdat)
  p <- pcf(simdat, spar=0.5, method="b")
  plot(p, main="pair correlation function for simdat")
  # indicates inhibition at distances r < 0.3

  # multitype point pattern
  data(betacells)
  p <- pcf(alltypes(betacells, "K"), spar=0.5, method="b")
  plot(p)
  # short range inhibition between all types
  # strong inhibition between "on" and "off"
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
