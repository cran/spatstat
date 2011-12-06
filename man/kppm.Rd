\name{kppm}
\alias{kppm}
\title{Fit Cluster or Cox Point Process Model}
\description{
  Fit a homogeneous or inhomogeneous cluster process or
  Cox point process model to a point pattern.
}
\usage{
  kppm(X, trend = ~1, clusters = "Thomas", covariates = NULL, ...,
       statistic="K", statargs=list())
}
\arguments{
  \item{X}{
    Point pattern (object of class \code{"ppp"}) to which the model
    should be fitted.
  }
  \item{trend}{
    An \R formula, with no left hand side,
    specifying the form of the log intensity.
  }
  \item{clusters}{
    Character string determining the cluster model.
    Partially matched.
    Options are \code{"Thomas"}, \code{"MatClust"},
    \code{"Cauchy"}, \code{"VarGamma"} and \code{"LGCP"}.
  }
  \item{covariates}{
    The values of any spatial covariates (other than the Cartesian
    coordinates) required by the model.
    A named list of pixel images, functions, windows or numeric constants.
  }
  \item{\dots}{
    Arguments passed to 
    \code{\link{thomas.estK}} or
    \code{\link{thomas.estpcf}} or
    \code{\link{matclust.estK}} or
    \code{\link{matclust.estpcf}} or
    \code{\link{lgcp.estK}} or
    \code{\link{lgcp.estpcf}} or
    \code{\link{cauchy.estK}} or
    \code{\link{cauchy.estpcf}} or
    \code{\link{vargamma.estK}} or
    \code{\link{vargamma.estpcf}}
    controlling the minimum contrast fitting algorithm.
  }
  \item{statistic}{
    The choice of summary statistic: either \code{"K"} or \code{"pcf"}.
  }
  \item{statargs}{
    Optional list of arguments to be used when calculating
    the summary statistic. See Details.
  }
}
\details{
  This function fits a Cox point process model to the
  point pattern dataset \code{X}. Cox models are suitable for
  spatially clustered point patterns.

  The model may be either a \emph{Poisson cluster process}
  or a \emph{Cox process}.
  The type of model is determined by the argument \code{clusters}.
  Currently the options 
  are \code{clusters="Thomas"} for the Thomas process,
  \code{clusters="MatClust"} for the Matern cluster process,
  \code{clusters="Cauchy"} for the Neyman-Scott cluster process
  with Cauchy kernel,
  \code{clusters="VarGamma"} for the Neyman-Scott cluster process
  with Variance Gamma kernel,
  and \code{clusters="LGCP"} for the log-Gaussian Cox process.
  The first four models are Poisson cluster processes.
  
  If the trend is constant (\code{~1})
  then the model is \emph{homogeneous}.
  The empirical \eqn{K}-function of the data is computed
  using \code{\link{Kest}},
  and the parameters of the cluster model are estimated by
  the method of minimum contrast (matching the theoretical
  \eqn{K}-function of the model to the empirical \eqn{K}-function
  of the data, as explained in \code{\link{mincontrast}}).

  Otherwise, the model is \emph{inhomogeneous}. 
  The algorithm first estimates the intensity function
  of the point process, by fitting a Poisson process with log intensity
  of the form specified by the formula \code{trend}.
  Then the inhomogeneous \eqn{K} function is estimated
  by \code{\link{Kinhom}} using this fitted intensity.
  Finally the parameters of the cluster model
  are estimated by the method of minimum contrast using the
  inhomogeneous \eqn{K} function. This two-step estimation
  procedure is due to Waagepetersen (2007).
  
  If \code{statistic="pcf"} then instead of using the
  \eqn{K}-function, the algorithm will use
  the pair correlation function \code{\link{pcf}} for homogeneous
  models and the inhomogeneous pair correlation function
  \code{\link{pcfinhom}} for inhomogeneous models.
  In this case, the smoothing parameters of the pair correlation
  can be controlled using the argument \code{statargs},
  as shown in the Examples.
}
\value{
  An object of class \code{"kppm"} representing the fitted model.
  There are methods for printing, plotting, predicting, simulating
  and updating objects of this class.
}
\seealso{
  methods for \code{kppm} objects:
  \code{\link{plot.kppm}},
  \code{\link{predict.kppm}},
  \code{\link{simulate.kppm}},
  \code{\link{update.kppm}},
  \code{\link{vcov.kppm}},
  \code{\link[spatstat:methods.kppm]{methods.kppm}},
  \code{\link{Kmodel.kppm}},
  \code{\link{pcfmodel.kppm}}.

  Fitting algorithms:
  \code{\link{thomas.estK}},
  \code{\link{matclust.estK}},
  \code{\link{lgcp.estK}},
  \code{\link{cauchy.estK}},
  \code{\link{vargamma.estK}},
  \code{\link{thomas.estpcf}},
  \code{\link{matclust.estpcf}},
  \code{\link{lgcp.estpcf}},
  \code{\link{cauchy.estpcf}},
  \code{\link{vargamma.estpcf}},
  \code{\link{mincontrast}}.

  Summary statistics:
  \code{\link{Kest}},
  \code{\link{Kinhom}},
  \code{\link{pcf}},
  \code{\link{pcfinhom}}.

  See also \code{\link{ppm}}
}
\references{
  Jalilian, A., Guan, Y. and Waagepetersen, R. (2011)
  Decomposition of variance for spatial Cox processes.
  Manuscript submitted for publication.

  Waagepetersen, R. (2007)
  An estimating function approach to inference for
  inhomogeneous Neyman-Scott processes.
  \emph{Biometrics} \bold{63}, 252--258.
}
\examples{
  data(redwood)
  kppm(redwood, ~1, "Thomas")
  kppm(redwood, ~x, "MatClust") 
  kppm(redwood, ~x, "MatClust", statistic="pcf", statargs=list(stoyan=0.2)) 
  kppm(redwood, ~1, "LGCP", statistic="pcf")
  kppm(redwood, ~x, cluster="Cauchy", statistic="K")
  kppm(redwood, cluster="VarGamma", nu.ker = 0.5, statistic="pcf")
  if(require(RandomFields)) {
     kppm(redwood, ~x, "LGCP", statistic="pcf",
           covmodel=list(model="matern", nu=0.3))
  }
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/},
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
  with contributions from Abdollah Jalilian and Rasmus Waagepetersen.
}
\keyword{spatial}
\keyword{models}

