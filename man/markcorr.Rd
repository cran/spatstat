\name{markcorr}
\alias{markcorr}
\title{
  Mark Correlation Function
}
\description{
  Estimate the marked correlation function
  of a marked point pattern.
}
\usage{
markcorr(X, f = function(m1,m2) { m1 * m2 }, r=NULL,
         slow=FALSE, correction=c("border", "isotropic", "Ripley", "translate"),
         method="density", \dots)
}
\arguments{
  \item{X}{The observed point pattern.
    An object of class \code{"ppp"} or something acceptable to
    \code{\link{as.ppp}}. 
  }
  \item{f}{Function \eqn{f} used in the definition of the
    mark correlation function.
  }
  \item{r}{numeric vector. The values of the argument \eqn{r}
    at which the mark correlation function 
    \eqn{\rho_f(r)}{rho_f(r)}
    should be evaluated.
    There is a sensible default.
  }
  \item{correction}{
    A character vector containing any selection of the
    options \code{"border"}, 
    \code{"isotropic"}, \code{"Ripley"} or \code{"translate"}.
    It specifies the edge correction(s) to be applied.
  }
  \item{method}{
    A character vector indicating the user's choice of
    density estimation technique to be used. Options are
    \code{"density"}, 
    \code{"loess"},
    \code{"sm"} and \code{"smrep"}.
  }
  
  \item{slow}{
    Logical vector indicating whether to use exact analytic geometry
    to calculate 
    the edge correction weights for the translation correction
    in the case where the window of observation of \code{X}
    is polygonal. These calculations are extremely slow and should
    be avoided unless you are a developer investigating the
    package.
  }
  \item{\dots}{
    Arguments passed to the density estimation routine
    (\code{\link{density}}, \code{\link{loess}} or \code{sm.density})
    selected by \code{method}.
  }
}
\value{
  An object of class \code{"fv"} (see \code{\link{fv.object}}).
  
  Essentially a data frame containing numeric columns 
  \item{r}{the values of the argument \eqn{r} 
    at which the mark correlation function \eqn{\rho_f(r)}{rho_f(r)}
    has been  estimated
  }
  \item{theo}{the theoretical value of \eqn{\rho_f(r)}{rho_f(r)}
    when the marks attached to different points are independent,
    namely 1
  }
  together with a column or columns named 
  \code{"border"}, 
  \code{"iso"} and/or \code{"trans"},
  according to the selected edge corrections. These columns contain
  estimates of the function \eqn{\rho_f(r)}{rho_f(r)}
  obtained by the edge corrections named.
}
\details{
  The mark correlation function \eqn{\rho_f(r)}{rho_f(r)}
  of a marked point process \eqn{X}
  is a measure of the dependence between the marks of two 
  points of the process a distance \eqn{r} apart.
  It is informally defined as
  \deqn{
    \rho_f(r) = \frac{E[f(M_1,M_2)]}{E[f(M,M')]}
  }{
    rho_f(r) = E[f(M1,M2)]/E[f(M,M')]
  }
  where \eqn{E[ ]} denotes expectation and \eqn{M_1,M_2}{M1,M2}
  are the marks attached to two points of the process
  separated by a distance \eqn{r}, while \eqn{M,M'} are
  independent realisations of the marginal distribution of marks.

  Here \eqn{f} is any function
  \eqn{f(m_1,m_2)}{f(m1,m2)}
  with two arguments which are possible marks of the pattern,
  and which returns a nonnegative real value.
  Common choices of \eqn{f} are:
  for continuous real-valued marks,
  \deqn{f(m_1,m_2) = m_1 m_2}{f(m1,m2)= m1 * m2}
  for discrete marks (multitype point patterns),
  \deqn{f(m_1,m_2) = 1(m_1 = m_2)}{f(m1,m2)= (m1 == m2)}
  and for marks taking values in \eqn{[0,2\pi)}{[0,2 * pi)},
  \deqn{f(m_1,m_2) = \sin(m_1 - m_2)}{f(m1,m2) = sin(m1-m2)}.
  
  Note that \eqn{\rho_f(r)}{rho_f(r)} is not a ``correlation''
  in the usual statistical sense. It can take any 
  nonnegative real value. The value 1 suggests ``lack of correlation'':
  if the marks attached to the points of \code{X} are independent
  and identically distributed, then
  \eqn{\rho_f(r) \equiv 1}{rho_f(r) =  1}.
  The interpretation of values larger or smaller than 1 depends
  on the choice of function \eqn{f}.

  The argument \code{X} must be a point pattern (object of class
  \code{"ppp"}) or any data that are acceptable to \code{\link{as.ppp}}.
  It must be a marked point pattern.

  The argument \code{f}
  must be a function, accepting two arguments \code{m1}
  and \code{m2} which are vectors of equal length containing mark
  values (of the same type as the marks of \code{X}).
  It must return a vector of numeric
  values of the same length as \code{m1} and \code{m2}.
  The values must be non-negative.

  The argument \code{r} is the vector of values for the
  distance \eqn{r} at which \eqn{\rho_f(r)}{rho_f(r)} is estimated.

  This algorithm assumes that \code{X} can be treated
  as a realisation of a stationary (spatially homogeneous) 
  random spatial point process in the plane, observed through
  a bounded window.
  The window (which is specified in \code{X} as \code{X$window})
  may have arbitrary shape.

  Biases due to edge effects are
  treated in the same manner as in \code{\link{Kest}}.
  The edge corrections implemented here are
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
  homogeneous). 

  The numerator and denominator of the mark correlation function
  (in the expression above) are estimated using density estimation
  techniques. The user can choose between
  \describe{
    \item{\code{"density"}}{
      which uses the standard kernel
      density estimation routine \code{\link{density}}, and
      works only for evenly-spaced \code{r} values;
    }
    \item{\code{"loess"}}{
      which uses the function \code{loess} in the
      package \code{modreg};
    }
    \item{\code{"sm"}}{
      which uses the function \code{sm.density} in the
      package \code{sm} and is extremely slow;
    }
    \item{\code{"smrep"}}{
      which uses the function \code{sm.density} in the
      package \code{sm} and is relatively fast.
    }
  }
}
\references{
  Stoyan, D. and Stoyan, H. (1994)
  Fractals, random shapes and point fields:
  methods of geometrical statistics.
  John Wiley and Sons.
}
\seealso{
 \code{\link{Kest}},
 \code{\link{Kmulti}}
}
\examples{
    # CONTINUOUS-VALUED MARKS:
    # (1) Longleaf Pine data
    # marks represent tree diameter
    data(longleaf)
    # Subset of this large pattern
    swcorner <- owin(c(0,100),c(0,100))
    sub <- longleaf[ , swcorner]
    # mark correlation function
    mc <- markcorr(sub)
    plot(mc)

    # (2) simulated data with independent marks
    X <- rpoispp(100)
    X <- X \%mark\% runif(X$n)
    Xc <- markcorr(X)
    plot(Xc)
    
    # MULTITYPE DATA:
    # Hughes' amacrine data
    # Cells marked as 'on'/'off'
    data(amacrine)
    # (3) Kernel density estimate with Epanecnikov kernel
    # (as proposed by Stoyan & Stoyan)
    M <- markcorr(amacrine, function(m1,m2) {m1==m2},
                  correction="translate", method="density",
                  kernel="epanechnikov")
    plot(M)
    # Note: kernel="epanechnikov" comes from help(density)

    # (4) Same again with explicit control over bandwidth
    M <- markcorr(amacrine, function(m1,m2) {m1==m2},
                  correction="translate", method="density",
                  kernel="epanechnikov", bw=0.02)
    # see help(density) for correct interpretation of 'bw'

   \testonly{
    data(ganglia)
    niets <- markcorr(ganglia, function(m1,m2){m1 == m2}, method="loess")
    niets <- markcorr(X, correction="isotropic", method="smrep")
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


