\name{mpl}
\alias{mpl}
\title{
  Fit Point Process Model by Maximum Pseudolikelihood
}
\description{
  Fits a point process model to an observed point pattern
  by the method of maximum pseudolikelihood.
}
\usage{
 mpl(Q, trend=~1, interaction=NULL, data, correction="border", rbord=0, use.gam=FALSE)
}
\arguments{
  \item{Q}{
    A data point pattern (of class \code{"ppp"})
    to which the model will be fitted,
    or preferably a quadrature scheme (of class \code{"quad"})
    containing this pattern.
  }
  \item{trend}{
  An \R formula object specifying the spatial trend to be fitted. 
  The default formula, \code{~1}, indicates the model is stationary
  and no trend is to be fitted. 
  }
  \item{interaction}{
    An object of class \code{"interact"}
    describing the point process interaction
    structure, or \code{NULL} indicating that a Poisson process (stationary
    or nonstationary) should be fitted.
  }
  \item{data}{
    An optional data frame of spatial covariates (evaluated at the
    locations given in the quadrature scheme \code{Q}).
  }
  \item{correction}{
    The name of the edge correction to be used. The default 
    is \code{"border"} indicating the border correction.
    Other possibilities may include \code{"Ripley"}, \code{"isotropic"},
    \code{"translate"} and \code{"none"}, depending on the 
    \code{interaction}.
  }
  \item{rbord}{
    If \code{correction = "border"}
    this argument specifies the distance by which
    the window should be eroded for the border correction.
  }
  \item{use.gam}{
    Logical flag; if \code{TRUE} then computations are performed
    using \code{gam} instead of \code{\link{glm}}.
  }
}
\value{
  An object of class \code{"ppm"} describing a fitted point process
  model.
 
  The fitted parameters can be obtained just by printing this object.
  Fitted spatial trends can be extracted using the \code{predict} method
  for this object (see \code{\link{predict.ppm}}).

  See \code{\link{ppm.object}} for details of the format of this object.
}
\details{
  This function fits a point process model
  to an observed point pattern by the method of maximum pseudolikelihood
  (Besag, 1975).
  The model may include
  spatial trend, interpoint interaction, and dependence on covariates.
 
  The algorithm is an implementation of the method of
  Baddeley and Turner (2000), which approximates the pseudolikelihood
  by a special type of quadrature sum generalising the Berman-Turner (1992) 
  approximation. 
 
  The argument \code{Q} should be either a point pattern
  or a quadrature scheme. If it is a point pattern, it is converted
  into a quadrature scheme.

  A quadrature scheme is an object of class \code{"quad"}
  (see \code{\link{quad.object}})
  which specifies both the data point pattern and the dummy points
  for the quadrature scheme, as well as the quadrature weights
  associated with these points.
  If \code{Q} is simply a point pattern
  (of class \code{"ppp"}, see \code{\link{ppp.object}})
  then it is interpreted as specifying the
  data points only; a set of dummy points specified
  by \code{\link{default.dummy}()} is added, and the default weighting rule is
  invoked to compute the quadrature weights.
 
  The usage of \code{mpl()} is closely analogous to the Splus/\R functions
  \code{\link{glm}} and \code{gam}.
  The analogy is:
  \tabular{ll}{
  \bold{glm} \tab \bold{mpl} \cr
  \code{formula} \tab \code{trend} \cr
  \code{family} \tab \code{interaction}
  }
  The point process model to be fitted is specified by the 
  arguments \code{trend} and \code{interaction}
  which are respectively analogous to
  the \code{formula} and \code{family} arguments of glm(). 
 
  Systematic effects (spatial trend and/or dependence on 
  spatial covariates) are specified by the argument
  \code{trend}. This is an Splus/\R formula object, which may be expressed
  in terms of the Cartesian coordinates \code{x}, \code{y},
  the marks \code{marks},
  or the variables in the data frame \code{data} (if supplied), or both.
  It specifies the \bold{logarithm} of the first order potential
  of the process.
  The formula should not use the names \code{Y}, \code{V}, \code{W}, or
  \code{SUBSET},
  which are reserved for internal use.
  If \code{trend} is absent or equal to the default, \code{~1}, then
  the model to be fitted is stationary (or at least, its first order 
  potential is constant). 
 
  Stochastic interactions between random points of the point process
  are defined by the argument \code{interaction}. This is an object of
  class \code{"interact"} which is initialised in a very similar way to the
  usage of family objects in \code{\link{glm}} and \code{gam}.
  See the examples below.
 
  If \code{interaction} is missing or \code{NULL}, then the model to be fitted
  has no interpoint interactions, that is, it is a Poisson process
  (stationary or nonstationary according to \code{trend}). In this case
  the method of maximum pseudolikelihood
  coincides with maximum likelihood.
 
  The argument \code{data}, if supplied, must be a data frame with
  as many rows as there are points in \code{Q}.
  The \eqn{i}th row of \code{data} should contain the values of
  spatial variables which have been observed
  at the \eqn{i}th point of \code{Q}. In this case \code{Q} must be a
  quadrature scheme, not merely a point pattern.
  Thus, it is not sufficient to have observed
  a spatial variable only at the points of the data point pattern; 
  the variable must also have been observed at certain other 
  locations in the window.
  The variable names \code{x}, \code{y} and \code{marks}
  are reserved for the Cartesian 
  coordinates and the mark values,
  and these should not be used for variables in \code{data}.
 
  The argument \code{correction} is the name of an edge correction method.
  The default \code{correction="border"} specifies the border correction,
  in which the quadrature window (the domain of integration of the 
  pseudolikelihood) is obtained by trimming off a margin of width \code{rbord}
  from the observation window of the data pattern.
  Not all edge corrections are implemented (or implementable)
  for arbitrary windows.
  Other options depend on the argument \code{interaction}, but these generally 
  include \code{"periodic"} (the periodic or toroidal edge correction
  in which opposite edges of a rectangular window are identified)
  and \code{"translate"}
  (the translation correction, see Baddeley 1998 and Baddeley and Turner
  2000).
  For pairwise interaction there is also Ripley's isotropic correction,
  identified by \code{"isotropic"} or \code{"Ripley"}.
 
  The fitted point process model returned by this function can be printed 
  (by the print method \code{\link{print.ppm}})
  to inspect the fitted parameter values.
  If a nonparametric spatial trend was fitted, this can be extracted using
  the predict method \code{\link{predict.ppm}}.

  This algorithm approximates the log pseudolikelihood
  by a sum over a finite set of quadrature points.
  Finer quadrature schemes (i.e. those with more
  quadrature points) generally yield a better approximation, at the
  expense of higher computational load. 
  Complete control over the quadrature scheme is possible.
  See \code{\link{quadscheme}} for an overview.

  Note that the method of maximum pseudolikelihood is
  believed to be inefficient and biased for point processes with strong
  interpoint interactions. In such cases, it is advisable to use
  iterative maximum likelihood
  methods such as Monte Carlo Maximum Likelihood (Geyer, 1999) 
  provided the appropriate simulation algorithm exists. 
  The maximum pseudolikelihood parameter estimate often serves
  as a good initial starting point for these iterative methods. 
  Maximum pseudolikelihood may also be used profitably for
  model selection in the initial phases of modelling.
}
\references{
  Baddeley, A. and Turner, R.
  Practical maximum pseudolikelihood for spatial point patterns.
  \emph{Australian and New Zealand Journal of Statistics}
  \bold{42} (2000) 283--322.
 
  Berman, M. and Turner, T.R. 
  Approximating point process likelihoods with GLIM.
  \emph{Applied Statistics} \bold{41} (1992) 31--38.
 
  Besag, J.
  Statistical analysis of non-lattice data.
  \emph{The Statistician} \bold{24} (1975) 179-195.
 
  Diggle, P.J., Fiksel, T., Grabarnik, P., Ogata, Y., Stoyan, D. and
  Tanemura, M.
  On parameter estimation for pairwise interaction processes.
  \emph{International Statistical Review} \bold{62} (1994) 99-117.
 
  Jensen, J.L. and Moeller, M.
  Pseudolikelihood for exponential family models of spatial point processes.
  \emph{Annals of Applied Probability} \bold{1} (1991) 445--461.
 
  Jensen, J.L. and Kuensch, H.R. 
  On asymptotic normality of pseudo likelihood
  estimates for pairwise interaction processes,
  \emph{Annals of the Institute of Statistical Mathematics}
  \bold{46} (1994) 475-486.
}

\seealso{
  \code{\link{ppp}}, 
  \code{\link{quadscheme}}, 
  \code{\link{ppm.object}},
  \code{\link{Poisson}},
  \code{\link{Strauss}},
  \code{\link{StraussHard}},
  \code{\link{Softcore}},
  \code{\link{Pairwise}},
  \code{\link{PairPiece}},
  \code{\link{Geyer}},
  \code{\link{Saturated}},
  \code{\link{OrdThresh}},
  \code{\link{Ord}}
}
\section{Warnings}{
See the comments above about the possible inefficiency
and bias of the maximum pseudolikelihood estimator.
 
The accuracy of the Berman-Turner-Baddeley approximation to
the pseudolikelihood depends on the number of dummy points used
in the quadrature scheme. The number of dummy points should 
at least equal the number of data points.
 
The parameter values of the fitted model
do not necessarily determine a valid point process.
Some of the point process models are only defined when the parameter
values lie in a certain subset. For example the Strauss process only 
exists when the interaction parameter \eqn{\gamma}{gamma}
is less than or equal to \eqn{1},
corresponding to a value of \code{mpl()$theta[2]}
less than or equal to \code{0}.
The current version of \code{mpl} maximises the pseudolikelihood
without constraining the parameters, and does not apply any checks for
sanity after fitting the model. 
 
The \code{trend} formula should not use the names
\code{Y}, \code{V}, \code{W}, or \code{SUBSET},
which are reserved
for internal use. The data frame \code{data} should have as many rows
as there are points in \code{Q}. It should not contain
variables called \code{x}, \code{y} or \code{marks}
as these names are reserved for the Cartesian coordinates
and the marks.
 
If the model formula involves one of the functions
\code{poly()}, \code{bs()}
or \code{ns()}
(e.g. applied to spatial coordinates \code{x} and \code{y}),
the fitted coefficients can be misleading.
The resulting fit is not to the raw spatial variates
(\code{x}, \code{x^2}, \code{x*y}, etc.) 
but to a transformation of these variates.  The transformation is implemented
by \code{poly()} in order to achieve better numerical stability.
However the
resulting coefficients are appropriate for use with the transformed
variates, not with the raw variates.  
This affects the interpretation of the constant
term in the fitted model, \code{logbeta}. 
Conventionally, \eqn{\beta}{beta} is the background intensity, i.e. the  
value taken by the conditional intensity function when all predictors
(including spatial or ``trend'' predictors) are set equal to \eqn{0}.
However the coefficient actually produced is the value that the
log conditional intensity takes when all the predictors, 
including the \emph{transformed}
spatial predictors, are set equal to \code{0}, which is not the same thing.

Worse still, the result of \code{\link{predict.ppm}} can be
completely wrong if the trend formula contains one of the
functions \code{poly()}, \code{bs()}
or \code{ns()}. This is a weakness of the underlying
function \code{\link{predict.glm}}. 

If you wish to fit a polynomial trend, 
we offer an alternative to \code{\link{poly}()},
namely \code{polynom()}, which avoids the
difficulty induced by transformations.  It is completely analogous
to \code{poly} except that it does not orthonormalise.
The resulting coefficient estimates then have
their natural interpretation and can be predicted correctly. 
Numerical stability may be compromised.

Values of the maximised pseudolikelihood are not comparable
if they have been obtained with different values of \code{rbord}.
}
\examples{
 \dontrun{library(spatstat)}
 data(nztrees)
 mpl(nztrees)
 # fit the stationary Poisson process
 # to point pattern 'nztrees'

 Q <- quadscheme(nztrees) 
 mpl(Q) 
 # equivalent.

 mpl(nztrees, ~ x)
 # fit the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(a + bx)
 # where x,y are the Cartesian coordinates
 # and a,b are parameters to be estimated

 mpl(nztrees, ~ polynom(x,2))
 # fit the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(a + bx + cx^2)

 library(splines)
 mpl(nztrees, ~ bs(x,df=3))
 #       WARNING: do not use predict.ppm() on this result
 # Fits the nonstationary Poisson process 
 # with intensity function lambda(x,y) = exp(B(x))
 # where B is a B-spline with df = 3
 
 mpl(nztrees, ~1, Strauss(r=10), rbord=10)
 # Fit the stationary Strauss process with interaction range r=10
 # using the border method with margin rbord=10
  
 mpl(nztrees, ~ x, Strauss(13), correction="periodic")
 # Fit the nonstationary Strauss process with interaction range r=13
 # and exp(first order potential) =  activity = beta(x,y) = exp(a+bx)
 # using the periodic correction.

 \dontrun{
 data(soilsurvey)
 mpl(soilsurvey, ~ bs(pH,3), Strauss(0.1), rbord=0.1, data=soilchem)
 #       WARNING: do not use predict.ppm() on this result
 # Fit the nonstationary Strauss process 
 # with intensity modelled as a third order spline function of the 
 # spatial variable "pH" in data frame 'soilchem'
 }

 ## MULTITYPE POINT PROCESSES ### 
 data(lansing)
 # Multitype point pattern --- trees marked by species
 mpl(lansing, ~ marks, Poisson())
 # fit stationary marked Poisson process
 # with different intensity for each species

\dontrun{
 mpl(lansing, ~ marks * polynom(x,y,3), Poisson())
}
 # fit nonstationary marked Poisson process
 # with different log-cubic trend for each species
\testonly{
 # equivalent functionality - smaller dataset
 data(ganglia)
 mpl(ganglia, ~ marks * polynom(x,y,2), Poisson())
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
