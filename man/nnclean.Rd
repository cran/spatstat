\name{nnclean}
\alias{nnclean}
\alias{nnclean.ppp}
\alias{nnclean.pp3}
\title{
  Nearest Neighbour Clutter Removal
}
\description{
  Detect features in a 2D or 3D spatial point pattern
  using nearest neighbour clutter removal.
}
\usage{
  nnclean(X, k, ...)
  \method{nnclean}{ppp}(X, k, ..., probs = FALSE,
                     edge.correct = FALSE, wrap = 0.1,
                     convergence = 0.001, plothist = FALSE,
                     verbose = TRUE, maxit = 50)
  \method{nnclean}{pp3}(X, k, ..., 
                     convergence = 0.001, plothist = FALSE,
                     verbose = TRUE, maxit = 50)
}
\arguments{
  \item{X}{
    A two-dimensional spatial point pattern (object of class
    \code{"ppp"}) or a three-dimensional point pattern
    (object of class \code{"pp3"}).
  }
  \item{k}{
    Degree of neighbour: \code{k=1} means nearest neighbour,
    \code{k=2} means second nearest, etc.
  }
  \item{\dots}{
    Ignored.
  }
  \item{probs}{
    Logical flag determining the type of value returned.
    See Details.
  }
  \item{edge.correct}{
    Logical flag specifying whether periodic edge correction
    should be performed (only implemented in 2 dimensions).
  }
  \item{wrap}{
    Numeric value specifying the relative size of the margin
    in which data will be replicated for the 
    periodic edge correction (if \code{edge.correct=TRUE}).
    A fraction of window width and window height.
  }
  \item{convergence}{
    Tolerance threshold for testing convergence of EM algorithm.
  }
  \item{maxit}{
    Maximum number of iterations for EM algorithm.
  }
  \item{plothist}{
    Logical flag specifying whether to plot a diagnostic histogram
    of the nearest neighbour distances and the fitted distribution.
  }
  \item{verbose}{
    Logical flag specifying whether to print progress reports.
  }
}
\details{
  Byers and Raftery (1998) developed a technique for recognising
  features in a spatial point pattern in the presence of
  random clutter.

  For each point in the pattern, the distance to the
  \eqn{k}th nearest neighbour is computed. Then the E-M algorithm is
  used to fit a mixture distribution to the nearest neighbour distances.
  The mixture components represent the feature and the clutter. The
  mixture model can be used to classify each point as belong to one
  or other component.

  The function \code{nnclean} is generic, with methods for
  two-dimensional point patterns (class \code{"ppp"})
  and three-dimensional point patterns (class \code{"pp3"})
  currently implemented.

  For two-dimensional point patterns,
  \itemize{
    \item 
    in the default case \code{probs=FALSE}, the result is a
    binary classification. Each point is classified by maximum likelihood
    as belonging either to the feature or the clutter.
    The return value is a marked point pattern, where the marks are a factor
    with levels \code{"noise"} and \code{"feature"} indicating the
    classification of each point. 

    If \code{probs=TRUE}, the result is a fuzzy or probabilistic
    classification. The return value is a marked point pattern
    with marks giving the estimated probabilities
    that each point belongs to a feature.
  }

  For three-dimensional point patterns, the result is a 3D point
  pattern with two (additional) columns of marks, one containing the
  probabilities, and one containing the binary classification.
}
\value{
  An object of the same kind as \code{X},
  obtained by attaching marks to the points of \code{X}.

  If \code{probs=FALSE} (the default), the marks are a factor
  with levels \code{"noise"} and \code{"feature"} indicating the
  maximum likelihood classification of each point by the algorithm.

  If \code{probs=TRUE}, the marks are numeric,
  giving the estimated probabilities
  that each point belongs to a feature. 
}
\references{
  Byers, S. and Raftery, A.E. (1998)
  Nearest-neighbour clutter removal for estimating features
  in spatial point processes.
  \emph{Journal of the American Statistical Association}
  \bold{93}, 577--584.
}
\author{
  Original by Simon Byers and Adrian Raftery.
  Adapted for \pkg{spatstat} by Adrian Baddeley.
}
\seealso{
  \code{\link{nndist}}, 
  \code{\link{split.ppp}}, 
  \code{\link{cut.ppp}}
}
\examples{
  data(shapley)
  X <- nnclean(shapley, k=17)
  plot(X, chars=c(".", "+"), cols=1:2)
  Y <- split(X)
  plot(Y, chars="+", cex=0.5)
  Z <- nnclean(shapley, k=17, probs=TRUE)
  plot(cut(Z, breaks=3), chars=c(".", "+", "+"), cols=1:3)
}
\keyword{spatial}
\keyword{classif}
