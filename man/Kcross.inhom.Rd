\name{Kcross.inhom}
\alias{Kcross.inhom}
\title{
  Inhomogeneous Cross K Function
}
\description{
  For a multitype point pattern, 
  estimate the inhomogeneous version of the cross \eqn{K} function,
  which counts the expected number of points of type \eqn{j}
  within a given distance of a point of type \eqn{i},
  adjusted for spatially varying intensity.
}
\usage{
Kcross.inhom(X, i=1, j=2, lambdaI, lambdaJ, \dots,  r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") ,
lambdaIJ=NULL)
}
\arguments{
  \item{X}{The observed point pattern, 
    from which an estimate of the inhomogeneous cross type \eqn{K} function
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
  \item{lambdaI}{
    Values of the the estimated intensity of the sub-process of
    points of type \code{i}.
    Either a pixel image (object of class \code{"im"}),
    or a numeric vector containing the intensity values
    at each of the type \code{i} points in \code{X}. 
  }
  \item{lambdaJ}{
    Values of the the estimated intensity of the sub-process of
    points of type \code{j}.
    Either a pixel image (object of class \code{"im"}),
    or a numeric vector containing the intensity values
    at each of the type \code{j} points in \code{X}. 
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
    \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{\dots}{
    Ignored.
  }
  \item{lambdaIJ}{
    Optional. A matrix containing estimates of the
    product of the intensities \code{lambdaI} and \code{lambdaJ}
    for each pair of points of types \code{i} and \code{j} respectively.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).

  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the function \eqn{K_{ij}(r)}{Kij(r)} has been  estimated
  }
  \item{theo}{the theoretical value of  \eqn{K_{ij}(r)}{Kij(r)}
    for a marked Poisson process, namely \eqn{\pi r^2}{pi * r^2}
  }
  together with a column or columns named 
  \code{"border"}, \code{"bord.modif"},
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{K_{ij}(r)}{Kij(r)}
  obtained by the edge corrections named.
}
\details{
  This is a generalisation of the function \code{\link{Kcross}}
  to include an adjustment for spatially inhomogeneous intensity,
  in a manner similar to the function \code{\link{Kinhom}}.

  The inhomogeneous cross-type \eqn{K} function is described by
  Moller and Waagepetersen (2003, pages 48-49 and 51-53).
  
  Briefly, given a multitype point process, suppose the sub-process
  of points of type \eqn{j} has intensity function
  \eqn{\lambda_j(u)}{lambda[j](u)} at spatial locations \eqn{u}.
  Suppose we place a mass of \eqn{1/\lambda_j(\zeta)}{1/lambda[j](z)}
  at each point \eqn{\zeta}{z} of type \eqn{j}. Then the expected total
  mass per unit area is 1. The
  inhomogeneous ``cross-type'' \eqn{K} function 
  \eqn{K_{ij}^{\mbox{inhom}}(r)}{K[ij]inhom(r)} equals the expected
  total mass within a radius \eqn{r} of a point of the process
  of type \eqn{i}. 
  
  If the process of type \eqn{i} points
  were independent of the process of type \eqn{j} points,
  then \eqn{K_{ij}^{\mbox{inhom}}(r)}{K[ij]inhom(r)}
  would equal \eqn{\pi r^2}{pi * r^2}.
  Deviations between the empirical \eqn{K_{ij}}{Kij} curve
  and the theoretical curve \eqn{\pi r^2}{pi * r^2} 
  suggest dependence between the points of types \eqn{i} and \eqn{j}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern, and the mark vector
  \code{X$marks} must be a factor.

  The arguments \code{i} and \code{j} will be interpreted as
  levels of the factor \code{X$marks}. (Warning: this means that
  an integer value \code{i=3} will be interpreted as the 3rd smallest level,
  not the number 3). 

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
  Similarly \code{lambdaJ} should contain
  estimated values of the intensity of the sub-process of points of
  type \code{j}. It may be either a pixel image or a numeric vector.
  
  The optional argument \code{lambdaIJ} is for advanced use only.
  It is a matrix containing estimated
  values of the products of these two intensities for each pair of
  data points of types \code{i} and \code{j} respectively.
  
  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{K_{ij}(r)}{Kij(r)} should be evaluated. 
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
  The arguments \code{i} and \code{j} are interpreted as
  levels of the factor \code{X$marks}. Beware of the usual
  trap with factors: numerical values are not
  interpreted in the same way as character values. 
}
\seealso{
 \code{\link{Kcross}},
 \code{\link{Kinhom}},
 \code{\link{pcf}}
}
\examples{
    # Lansing Woods data
    data(lansing)
    lansing <- lansing[seq(1,lansing$n, by=10)]
    ma <- split(lansing)$maple
    wh <- split(lansing)$whiteoak

    # method (1): estimate intensities by nonparametric smoothing
    lambdaM <- density.ppp(ma, sigma=0.15)
    lambdaW <- density.ppp(wh, sigma=0.15)
    K <- Kcross.inhom(lansing, "whiteoak", "maple", lambdaW[wh], lambdaM[ma])
    K <- Kcross.inhom(lansing, "whiteoak", "maple", lambdaW,     lambdaM)

    # method (2): fit parametric intensity model
    fit <- ppm(lansing, ~marks * polynom(x,y,2))
    # evaluate fitted intensities at data points
    # (these are the intensities of the sub-processes of each type)
    inten <- fitted(fit, dataonly=TRUE)
    # split according to types of points
    lambda <- split(inten, lansing$marks)
    K <- Kcross.inhom(lansing, "whiteoak", "maple",
              lambda$whiteoak, lambda$maple)
    
    # synthetic example: type A points have intensity 50,
    #                    type B points have intensity 100 * x
    lamB <- as.im(function(x,y){50 + 100 * x}, owin())
    X <- superimpose(A=runifpoispp(50), B=rpoispp(lamB))
    XA <- split(X)$A
    XB <- split(X)$B
    K <- Kcross.inhom(X, "A", "B",
        lambdaI=rep(50, XA$n), lambdaJ=lamB[XB])
    K <- Kcross.inhom(X, "A", "B",
        lambdaI=as.im(50, X$window), lambdaJ=lamB)
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{nonparametric}
