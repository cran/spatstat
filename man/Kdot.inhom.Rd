\name{Kdot.inhom}
\alias{Kdot.inhom}
\title{
  Inhomogeneous Multitype K Dot Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the dot \eqn{K} function,
  which counts the expected number of points of any type
  within a given distance of a point of type \eqn{i},
  adjusted for spatially varying intensity.
}
\usage{
Kdot.inhom(X, i, lambdaI, lambdadot, \dots, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
         lambdaIdot=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{K} function
    \eqn{K_{i\bullet}(r)}{Ki.(r)} will be computed.
    It must be a multitype point pattern (a marked point pattern
    whose marks are a factor). See under Details.
  }
  \item{i}{Number or character string identifying the type (mark value)
    of the points in \code{X} from which distances are measured.
    Defaults to the first level of \code{marks(X)}.
  }
  \item{lambdaI}{
    Values of the the estimated intensity of the sub-process of
    points of type \code{i}.
    Either a pixel image (object of class \code{"im"}),
    or a numeric vector containing the type \code{i} intensity values
    at each of the type \code{i} points in \code{X}. 
  }
  \item{lambdadot}{
    Values of the the estimated intensity of the entire point process,
    Either a pixel image (object of class \code{"im"}),
    or a numeric vector containing the intensity values
    at each point in \code{X}. 
  }
  \item{\dots}{
    Ignored.
  }
  \item{r}{
      Optional. Numeric vector giving the values of the argument \eqn{r}
      at which the cross K function
      \eqn{K_{ij}(r)}{Kij(r)} should be evaluated.
      There is a sensible default.
      First-time users are strongly advised not to specify this argument.
      See below for important conditions on \eqn{r}.
  }
  \item{breaks}{
    Optional. An alternative to the argument \code{r}.
    Not normally invoked by the user. See the \bold{Details} section.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, \code{"bord.modif"},
    \code{"isotropic"}, \code{"Ripley"}, \code{"translate"},
    \code{"none"} or \code{"best"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{lambdaIdot}{
    Optional. A matrix containing estimates of the
    product of the intensities \code{lambdaI} and \code{lambdadot}
    for each pair of points, the first point of type \code{i} and
    the second of any type.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{K_{i\bullet}(r)}{Ki.(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{K_{i\bullet}(r)}{Ki.(r)}
    for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K_{i\bullet}(r)}{Ki.(r)}
  obtained by the edge corrections named.
}
\details{
  This is a generalisation of the function \code{\link{Kdot}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Kinhom}}.

  Briefly, given a multitype point process, consider the points without
  their types, and suppose this unmarked point process 
  has intensity function
  \eqn{\lambda(u)}{lambda(u)} at spatial locations \eqn{u}.
  Suppose we place a mass of \eqn{1/\lambda(\zeta)}{1/lambda(z)}
  at each point \eqn{\zeta}{z} of the process. Then the expected total
  mass per unit area is 1. The
  inhomogeneous ``dot-type'' \eqn{K} function 
  \eqn{K_{i\bullet}^{\mbox{inhom}}(r)}{K[i.]inhom(r)} equals the expected
  total mass within a radius \eqn{r} of a point of the process
  of type \eqn{i}, discounting this point itself.
  
  If the process of type \eqn{i} points
  were independent of the points of other types,
  then \eqn{K_{i\bullet}^{\mbox{inhom}}(r)}{K[i.]inhom(r)}
  would equal \eqn{\pi r^2}{pi * r^2}.
  Deviations between the empirical \eqn{K_{i\bullet}}{Ki.} curve
  and the theoretical curve \eqn{\pi r^2}{pi * r^2} 
  suggest dependence between the points of types \eqn{i} and \eqn{j} for
  \eqn{j\neq i}{j != i}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern, and the mark vector
  \code{X$marks} must be a factor.

  The argument \code{i} will be interpreted as a
  level of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the 3rd smallest level,
  not the number 3).
  If \code{i} is missing, it defaults to the first
  level of the marks factor, \code{i = levels(X$marks)[1]}.
  
  The argument \code{lambdaI} supplies the values
  of the intensity of the sub-process of points of type \code{i}.
  It may be either
  \describe{
    \item{a pixel image}{(object of class \code{"im"}) which
      gives the values of the type \code{i} intensity
      at all locations in the window containing \code{X};
    }
    \item{a numeric vector}{containing the values of the
      type \code{i} intensity evaluated only
      at the data points of type \code{i}. The length of this vector
      must equal the number of type \code{i} points in \code{X}.
    }
  }
  The argument \code{lambdadot} should contain
  estimated values of the intensity of the entire point process.
  It may be either a pixel image, or a numeric vector of length equal
  to the number of points in \code{X}.

  For advanced use only, the optional argument \code{lambdaIdot}
  is a matrix containing estimated
  values of the products of these two intensities for each pair of
  points, the first point of type \code{i} and the second of any type.
  
  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{i\bullet}(r)}{Ki.(r)} should be evaluated. 
  The values of \eqn{r} must be increasing nonnegative numbers
  and the maximum \eqn{r} value must exceed the radius of the
  largest disc contained in the window.

  The argument \code{correction} chooses the edge correction
  as explained e.g. in \code{\link{Kest}}.

  The pair correlation function can also be applied to the
  result of \code{Kcross.inhom}; see \code{\link{pcf}}.
}
\references{
  Moller, J. and Waagepetersen, R.
  Statistical Inference and Simulation for Spatial Point Processes
  Chapman and Hall/CRC
  Boca Raton, 2003.
}
\section{Warnings}{
  The argument \code{i} is interpreted as a
  level of the factor \code{X$marks}. Beware of the usual
  trap with factors: numerical values are not
  interpreted in the same way as character values. 
}
\seealso{
 \code{\link{Kdot}},
 \code{\link{Kinhom}},
 \code{\link{Kcross.inhom}},
 \code{\link{pcf}}
}
\examples{
    # Lansing Woods data
    data(lansing)
    lansing <- lansing[seq(1,lansing$n, by=10)]
    ma <- split(lansing)$maple
    lg <- unmark(lansing)

    # Estimate intensities by nonparametric smoothing
    lambdaM <- density.ppp(ma, sigma=0.15)
    lambdadot <- density.ppp(lg, sigma=0.15)
    K <- Kdot.inhom(lansing, "maple", lambdaI=lambdaM[ma],
                                      lambdadot=lambdadot[lg])
    K <- Kdot.inhom(lansing, "maple", lambdaI=lambdaM,
                                      lambdadot=lambdadot)


    # synthetic example: type A points have intensity 50,
    #                    type B points have intensity 50 + 100 * x
    lamB <- as.im(function(x,y){50 + 100 * x}, owin())
    lamdot <- as.im(function(x,y) { 100 + 100 * x}, owin())
    X <- superimpose(A=runifpoispp(50), B=rpoispp(lamB))
    XB <- split(X)$B
    K <- Kdot.inhom(X, "B",  lambdaI=lamB[XB], lambdadot=lamdot[X])
    K <- Kdot.inhom(X, "B",  lambdaI=lamB,     lambdadot=lamdot)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}

