\name{Kcross}
\alias{Kcross}
\title{
  Multitype K Function (Cross-type)
}
\description{
  For a multitype point pattern, 
  estimate the multitype \eqn{K} function
  which counts the expected number of points of type \eqn{j}
  within a given distance of a point of type \eqn{i}.
}
\synopsis{
Kcross(X, i=1, j=2, r=NULL, breaks=NULL, \dots)
}
\usage{
Kcross(X, i=1, j=2)
Kcross(X, i=1, j=2, r)
Kcross(X, i=1, j=2, breaks)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the cross type \eqn{K} function
    \eqn{K_{ij}(r)}{Kij(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
  }
  \item{j}{Number or character string identifying the type (mark value)
    of the points in \code{X} to which distances are measured.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the distribution function
    \eqn{K_{ij}(r)}{Kij(r)} should be evaluated.
    There is a sensible default.
    First-time users are strongly advised not to specify this argument.
    See below for important conditions on \eqn{r}.
  }
  \item{breaks}{An alternative to the argument \code{r}.
    Not normally invoked by the user. See the \bold{Details} section.
  }
}
\value{
  A data frame containing 3 numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{K_{ij}(r)}{Kij(r)} has been  estimated
  }
  \item{border}{the ``border correction''
    estimator of \eqn{K_{ij}(r)}{Kij(r)}
  }
  \item{theo}{the theoretical value of  \eqn{K_{ij}(r)}{Kij(r)}
    for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
  }
}
\details{
  This function \code{Kcross} and its companions
  \code{\link{Kdot}} and \code{\link{Kmulti}}
  are generalisations of the function \code{\link{Kest}}
  to multitype point patterns. 

  A multitype point pattern is a spatial pattern of
  points classified into a finite number of possible
  ``colours'' or ``types''. In the \code{spatstat} package,
  a multitype pattern is represented as a single 
  point pattern object in which the points carry marks,
  and the mark value attached to each point
  determines the type of that point.
  
  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern, and the mark vector
  \code{X$marks} must be a factor.
  The arguments \code{i} and \code{j} will be interpreted as
  levels of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the 3rd smallest level,
  not the number 3). 
  
  The ``cross-type'' (type \eqn{i} to type \eqn{j})
  \eqn{K} function 
  of a stationary multitype point process \eqn{X} is defined so that
  \eqn{\lambda_j K_{ij}(r)}{lambda[j] Kij(r)} equals the expected number of
  additional random points of type \eqn{j}
  within a distance \eqn{r} of a
  typical point of type \eqn{i} in the process \eqn{X}.
  Here \eqn{\lambda_j}{lambda[j]}
  is the intensity of the type \eqn{j} points,
  i.e. the expected number of points of type \eqn{j} per unit area.
  The function \eqn{K_{ij}}{Kij} is determined by the 
  second order moment properties of \eqn{X}.

  An estimate of \eqn{K_{ij}(r)}{Kij(r)}
  is a useful summary statistic in exploratory data analysis
  of a multitype point pattern.
  If the process of type \eqn{i} points
  were independent of the process of type \eqn{j} points,
  then \eqn{K_{ij}(r)}{Kij(r)} would equal \eqn{\pi r^2}{pi * r^2}.
  Deviations between the empirical \eqn{K_{ij}}{Kij} curve
  and the theoretical curve \eqn{\pi r^2}{pi * r^2} 
  may suggest dependence between the points of types \eqn{i} and \eqn{j}.

  This algorithm estimates the distribution function \eqn{K_{ij}(r)}{Kij(r)} 
  from the point pattern \code{X}. It assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.
  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kest}},
  using the border correction.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{ij}(r)}{Kij(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must exceed the radius of the
  largest disc contained in the window.

  The pair correlation function can also be applied to the
  result of \code{Kcross}; see \code{\link{pcf}}.
}
\references{
  Cressie, N.A.C. \emph{Statistics for spatial data}.
    John Wiley and Sons, 1991.

  Diggle, P.J. \emph{Statistical analysis of spatial point patterns}.
    Academic Press, 1983.

  Harkness, R.D and Isham, V. (1983)
  A bivariate spatial point pattern of ants' nests.
  \emph{Applied Statistics} \bold{32}, 293--303
 
  Lotwick, H. W. and Silverman, B. W. (1982).
  Methods for analysing spatial processes of several types of points.
  \emph{J. Royal Statist. Soc. Ser. B} \bold{44}, 406--413.

  Ripley, B.D. \emph{Statistical inference for spatial processes}.
  Cambridge University Press, 1988.

  Stoyan, D, Kendall, W.S. and Mecke, J.
  \emph{Stochastic geometry and its applications}.
  2nd edition. Springer Verlag, 1995.

}
\section{Warnings}{
  The arguments \code{i} and \code{j} are interpreted as
  levels of the factor \code{X$marks}. Beware of the usual
  trap with factors: numerical values are not
  interpreted in the same way as character values. See the first example.

  The reduced sample estimator of \eqn{K_{ij}}{Kij} is pointwise approximately 
  unbiased, but need not be a valid distribution function; it may 
  not be a nondecreasing function of \eqn{r}. Its range is always 
  within \eqn{[0,1]}.
}
\seealso{
 \code{\link{Kdot}},
 \code{\link{Kest}},
 \code{\link{Kmulti}},
 \code{\link{pcf}}
}
\examples{
    \dontrun{library(spatstat)}
    data(catWaessle)
     # cat retina data: types 0,1 represent off/on
    K01 <- Kcross(catWaessle, "0", "1") # note: "0" not 0
    \dontrun{plot(K01$r, K01$border,
                  xlab="r", ylab="K01(r)",
                  type="l", ylim=c(0,1))}
    K10 <- Kcross(catWaessle, "1", "0")
    # Poisson theoretical curve
    \dontrun{lines(K01$r, K01$theo, lty=2)}

    # synthetic example    
    pp <- runifpoispp(50)
    pp$marks <- sample(c("a","b"), pp$n, replace=TRUE)
    K <- Kcross(pp, "a", "b")
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
