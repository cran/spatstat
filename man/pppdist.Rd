\name{pppdist}
\alias{pppdist}
\title{Optimal Match Between Two Point Patterns}
\description{
  Given two point patterns, find the optimal match between them.
}
\usage{
  pppdist(X, Y, q = 1, precision = 7, show.rprimal = FALSE, belowone = TRUE, timelag = 0)
}
\arguments{
  \item{X,Y}{Two point patterns (objects of class \code{"ppp"}).}
  \item{q}{
    Exponent of the Wasserstein distance
    (or \code{Inf} for the Prohorov distance).
  }
  \item{precision}{
    Index controlling accuracy of algorithm.
    Interpoint distances will be rounded to the nearest multiple of
    \code{10^(-precision)}.
  }
  \item{show.rprimal}{
    Logical. Whether to display a plot showing the
    iterative solution of the restricted primal problem.
  }
  \item{belowone}{
    Logical. Experimental use only.
    Indicates whether to rescale the distances by a fudge factor.
  }
  \item{timelag}{
    Time lag, in seconds, between successive displays of the
    iterative solution of the restricted primal problem.
  }
}
\details{
  Finds the matching between the point patterns \code{X} and \code{Y}
  which minimises the sum of the distances between matched points
  (if \code{q=1}), the maximum distance between matched points
  (if \code{q=Inf}), and in general the \code{1/q}th power of the sum of
  the \code{q}th powers of the distances between matched points.
  If \eqn{q < 1} this is known as the Wasserstein distance,
  and if \eqn{q=Inf}{q = \infty} it is the Prohorov distance.

  For finite exponents \code{q}, there is a fast C algorithm,
  which will handle patterns of 100 points without difficulty,
  but should not be used with thousands of points.
  If \code{show.rprimal=TRUE}, slower interpreted code is used
  to demonstrate the algorithm.
  For \code{q=Inf}, even slower interpreted R code is used,
  and this works only for very small point patterns.
}
\value{
  An object of class \code{pppmatching} that represents the
  matching. There are methods for \code{plot}, \code{print} and \code{summary}
  for this class.
}
\author{
  Dominic Schuhmacher.
}
\seealso{}
\examples{
   X <- runifpoint(42)
   Y <- runifpoint(42)
   pppdist(X, Y)
   pppdist(X[1:10], Y[1:10], q=Inf)
}
\keyword{spatial}
\keyword{math}
