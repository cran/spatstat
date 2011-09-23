\name{default.rmhcontrol}
\alias{default.rmhcontrol}
\title{Set Default Control Parameters for Metropolis-Hastings Algorithm.}
\description{
  Given a fitted point process model, this command
  sets appropriate default values of the
  parameters controlling the iterative behaviour
  of the Metropolis-Hastings algorithm.
}
\usage{
   default.rmhcontrol(model,\dots)
}
\arguments{
  \item{model}{
    A fitted point process model (object of class \code{"ppm"})
  }
  \item{\dots}{
    Additional arguments passed to \code{\link{rmhcontrol}}.
  }
}
\value{
  An object of class \code{"rmhcontrol"}. See \code{\link{rmhcontrol}}.
}
\details{
  This function sets default values of 
  the parameters controlling the iterative behaviour
  of the Metropolis-Hastings
  simulation algorithm, using defaults
  that would be appropriate for the fitted point
  process model \code{model}.

  See \code{\link{rmhcontrol}} for the full list of control parameters.
  
  The defaults set by \code{default.rmhcontrol} are as follows.
  The expansion parameter \code{expand} is set to
  \code{\link{default.expand}(model)}.
  The periodic distance flag \code{periodic} is set to
  \code{TRUE} if \code{as.owin(model)} is a rectangle,
  and \code{FALSE} otherwise.
  All other parameters revert to their defaults given in
  \code{\link{rmhcontrol.default}}.
  
  Additional control parameters can be passed through the
  \code{\dots} argument.
}
\seealso{
  \code{\link{rmhcontrol}},
  \code{\link{ppm}},
  \code{\link{default.expand}}
}
\examples{
  data(cells)
  fit <- ppm(cells, ~1, Strauss(0.1))
  default.rmhcontrol(fit)
}
\author{Adrian Baddeley
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\keyword{spatial}
\keyword{datagen}
