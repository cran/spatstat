\name{plot.ppm}
\alias{plot.ppm}
\title{plot a Fitted Point Process Model}
\description{
  Given a fitted point process model obtained by \code{\link{ppm}},
  plot the spatial trend and conditional intensity of the model.
}
\usage{
 plot.ppm(x, ngrid = c(40,40), superimpose = TRUE, 
                  trend = TRUE, cif = TRUE, pause = TRUE,
                  how=c("persp","image", "contour"), \dots)
}
\arguments{
  \item{x}{
    A fitted point process model, typically obtained from
    the model-fitting algorithm \code{\link{ppm}}.
    An object of class \code{"ppm"}.
  }
  \item{ngrid}{
    the spatial trend and conditional intensity are to be evaluated and
    displayed on a rectangular grid (\code{ngrid[1]} by \code{ngrid[2]}
    or \code{ngrid} by \code{ngrid})
    of points.
  }
  \item{superimpose}{
    logical flag; if \code{TRUE}, the original data point pattern
    will be superimposed on the plots.
  }
  \item{trend}{
    logical flag; if \code{TRUE}, the spatial trend will be plotted.
  }
  \item{cif}{
    logical flag; if \code{TRUE}, the conditional intensity will be
    plotted.
  }
  \item{pause}{
    logical flag indicating whether to pause with a prompt 
    after each plot. Set \code{pause=FALSE} if plotting to a file.
  }
  \item{how}{
    character string or character vector indicating the style or styles of
    plots to be performed.
  }
  \item{\dots}{
    extra arguments to the plotting functions
    \code{\link{persp}}, \code{\link{image}} and \code{\link{contour}}.
  }
}
\value{
  none.
}
\details{
  This is the \code{plot} method for the class \code{"ppm"}
  (see \code{\link{ppm.object}} for details of this class).
  
  It invokes \code{\link{predict.ppm}} to compute the spatial 
  trend and conditional intensity of the fitted point process model.
  See \code{\link{predict.ppm}} for more explanation about spatial trend
  and conditional intensity.
 
  The spatial trend
  and/or the conditional intensity of the fitted spatial point process model
  \code{x} are computed at the points of a rectangular grid, and these
  values are plotted successively using
  \code{\link{persp}}, \code{\link{image}} and \code{\link{contour}}
  (or only a selection of these three, if \code{how} is given).
  Extra graphical parameters controlling the display may be passed
  directly via the arguments \code{...} or
  indirectly reset using \code{\link{spatstat.options}}.

  The default action is to create a rectangular 40 by 40 grid of points
  in the observation window of the data point pattern, and evaluate
  the spatial trend and conditional intensity at these locations.
 
  If the fitted model had no spatial trend, then the default is 
  to omit plotting this (flat) surface, unless \code{trend=TRUE}
  is set explicitly.
 
  If the fitted model was Poisson, so that there were no spatial interactions,
  then the conditional intensity and spatial trend are identical, and the
  default is to omit the conditional intensity, unless \code{cif=TRUE} is set
  explicitly.

  If the fitted model was a marked point process, then 
  predictions are made and plotted for each possible mark value in turn.

}
\seealso{
  \code{\link{ppm}},
  \code{\link{ppm.object}},
  \code{\link{predict.ppm}},
  \code{\link{print.ppm}},
  \code{\link{persp}},
  \code{\link{image}},
  \code{\link{contour}},
  \code{\link{plot}},
  \code{\link{spatstat.options}}
}
\section{Warnings}{
  See warnings in \code{\link{predict.ppm}}.
}
\examples{
 \dontrun{
 data(cells)
 Q <- quadscheme(cells)
 m <- ppm(Q, ~1, Strauss(0.05))
 plot(m)
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
