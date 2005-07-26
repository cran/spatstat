\name{diagnose.ppm}
\alias{diagnose.ppm}
\title{
  Diagnostic Plots for Fitted Point Process Model
}
\description{
  Given a point process model fitted to a point pattern,
  produce diagnostic plots based on residuals.
}
\synopsis{
  diagnose.ppm(object, \dots, type="raw", which="all", sigma=NULL, 
               rbord=reach(object), cumulative=TRUE,
               plot.it=TRUE, rv = NULL, compute.sd=TRUE,
               compute.cts=TRUE, typename)
}
\usage{
  diagnose.ppm(object, \dots, type="raw", which="all", sigma=NULL, 
               rbord=reach(object), cumulative=TRUE,
               plot.it=TRUE,
               plot.neg="image", plot.smooth="imagecontour", plot.sd,
               rv = NULL, plot.sd, spacing=0.1,
               srange=NULL, monochrome=FALSE)
}
\arguments{
  \item{object}{
    The fitted point process model (an object of class \code{"ppm"})
    for which diagnostics should be produced. This object
    is usually obtained from \code{\link{ppm}}.
  }
  \item{type}{
    String indicating the type of residuals or weights to be used.
    Current options are \code{"eem"}
    for the Stoyan-Grabarnik exponential energy weights,
    \code{"raw"} for the raw residuals,
    \code{"inverse"} for the inverse-lambda residuals,
    and \code{"pearson"} for the Pearson residuals.
    A partial match is adequate.
  }
  \item{which}{
    Character string or vector indicating the choice(s) of
    plots to be generated. Options are
    \code{"all"}, \code{"marks"}, \code{"smooth"},
    \code{"x"}, \code{"y"} and \code{"sum"}.
    Multiple choices may be given but must be matched exactly.
    See Details.
  }
  \item{sigma}{
    Bandwidth for kernel smoother in \code{"smooth"} option.
  }
  \item{\dots}{
    Extra arguments, controlling either the resolution of the smoothed image
    (passed to \code{ksmooth.ppp}) or the appearance of the plots
    (passed to \code{plot.default}).
  }
  \item{rbord}{
    Width of border to avoid edge effects.
    The diagnostic calculations
    will be confined to those points of the data pattern which are
    at least \code{rbord} units away from the edge of the window.
  }
  \item{cumulative}{
    Logical flag indicating whether the lurking variable plots
    for the \eqn{x} and \eqn{y} coordinates will be the plots of
    cumulative sums of marks (\code{cumulative=TRUE}) or the
    plots of marginal integrals of the smoothed residual field
    (\code{cumulative=FALSE}).
  }
  \item{plot.it}{
    Logical value indicating whether 
    plots should be shown. If \code{plot.it=FALSE}, 
    the computed diagnostic quantities are returned without plotting them.
  }
  \item{plot.neg}{
    One of \code{"discrete"} or \code{"image"}
    indicating how the density part of the residual measure should be plotted.
  }
  \item{plot.smooth}{
    One of \code{"image"}, \code{"persp"}, \code{"contour"}
    or \code{"imagecontour"}
    indicating how the smoothed residual field should be plotted.
  }
  \item{plot.sd}{
    Logical value indicating whether 
    error bounds should be added to the \code{"x"} and \code{"y"}
    plots. The default is \code{TRUE} for Poisson models and
    \code{FALSE} for non-Poisson models. See Details.
  }
  \item{rv}{
    Usually absent. 
    If this argument is present, the values of the residuals will not be
    calculated from the fitted model \code{object} but will instead
    be taken directly from this vector.
  }
  \item{spacing}{
    The spacing between plot panels (when a four-panel plot
    is generated) expressed as a fraction of the width of the
    window of the point pattern.
  }
  \item{srange}{
    Vector of length 2 that will be taken as giving the range of values
    of the smoothed residual field, when generating an image plot of this
    field. This is useful if you want to generate diagnostic plots
    for two different fitted models using the same colour map. 
   }
  \item{monochrome}{
    Flag indicating whether images should be displayed in
    greyscale (suitable for publication) or in colour (suitable
    for the screen). The default is to display in colour.
   }
}
\value{
  An object of class \code{"diagppm"} which contains
  the coordinates needed to reproduce the selected plots.
  This object can be plotted using \code{plot.diagppm}
  and printed using \code{print.diagppm}.
}
\details{
  This function generates several diagnostic plots for a
  fitted point process model.
  The plots display the residuals from the fitted model
  (Baddeley et al, 2005)
  or alternatively the `exponential energy marks' (Stoyan and Grabarnik, 1991).
  These plots can be used to
  assess goodness-of-fit, to identify outliers in the data,
  and to reveal departures from the fitted model.
  See also the companion function \code{\link{qqplot.ppm}}.

  The argument \code{object} must be a fitted point process model
  (object of class \code{"ppm"}) typically produced by the maximum
  pseudolikelihood fitting algorithm \code{\link{ppm}}).

  The argument \code{type} selects the type of residual or weight
  that will be computed. Current options are:

  \describe{
    \item{\code{"eem"}:}{
    exponential energy marks (Stoyan and Grabarnik, 1991) 
    computed by \code{\link{eem}}.
    These are positive weights attached to the data points
    (i.e. the points of the point pattern dataset
    to which the model was fitted).
    If the fitted model is correct, then the sum of these weights
    for all data points in a spatial region \eqn{B}
    has expected value equal to the
    area of \eqn{B}. See \code{\link{eem}} for further explanation.
  }
  \item{\code{"raw"}, \code{"inverse"} or \code{"pearson"}:}{
    point process residuals (Baddeley et al, 2005)
    computed by the function \code{\link{residuals.ppm}}.
    These are residuals attached both to the data points and to some
    other points in the window of observation (namely, to the dummy
    points of the quadrature scheme used to fit the model).
    If the fitted model is correct, then the sum of the
    residuals in a spatial region \eqn{B} has mean zero.
    The options are
    \itemize{
      \item
      \code{"raw"}: the raw residuals;
      \item
      \code{"inverse"}: the `inverse-lambda' residuals,
      a counterpart of the exponential energy weights;
      \item
      \code{"pearson"}: the Pearson residuals.
    }
    See \code{\link{residuals.ppm}} for further explanation.
  }
  }

  The argument \code{which} selects the type of plot that is
  produced. Options are:
  \describe{
    \item{\code{"marks"}:}{
      plot the residual measure.
      For the exponential energy weights (\code{type="eem"})
      this displays circles centred at the points of the data pattern,
      with radii proportional to the exponential energy weights.
      For the residuals (\code{type="raw"}, \code{type="inverse"}
      or \code{type="pearson"}) this again displays circles centred at
      the points of the data pattern with radii proportional to the
      (positive) residuals, while the plotting of the negative residuals
      depends on the argument \code{plot.neg}. If
      \code{plot.neg="image"} then the negative part of the residual
      measure, which is a density, is plotted as a colour image.
      If \code{plot.neg="discrete"} then the discretised negative
      residuals (obtained by approximately integrating the negative
      density using the quadrature scheme of the fitted model)
      are plotted as squares centred at the dummy points
      with side lengths proportional to the (negative) residuals.
    }
    \item{\code{"smooth"}:}{
      plot a kernel-smoothed version of the residual measure.
      Each data or dummy point is taken to have a `mass' equal to its
      residual or exponential energy weight.
      (Note that residuals can be negative).
      This point mass is then replaced by
      a bivariate isotropic Gaussian density
      with standard deviation \code{sigma}.
      The value of the smoothed residual field at
      any point in the window is the sum of these weighted densities.
      If the fitted model is correct, this smoothed field
      should be flat, and its height should be close to 0
      (for the residuals) or 1 (for the exponential energy weights).
      The field is plotted either as an image, contour plot or
      perspective view of a surface, according to the
      argument \code{plot.smooth}.
      The range of values of the smoothed field is printed
      if the option \code{which="sum"} is also selected.
    }
    \item{\code{"x"}:}{
      produce a `lurking variable' plot for the \eqn{x} coordinate.
      This is a plot of \eqn{h(x)} against \eqn{x} (solid lines)
      and of \eqn{E(h(x))} against \eqn{x} (dashed lines),
      where \eqn{h(x)} is defined below, and \eqn{E(h(x))} denotes the
      expectation of \eqn{h(x)} assuming the fitted model is true.
      \itemize{
      \item
        if \code{cumulative=TRUE} then \eqn{h(x)} is the cumulative sum of
	the weights or residuals for all points
	which have \eqn{X} coordinate less than or equal to \eqn{x}.
	For the residuals \eqn{E(h(x)) = 0},
	and for the exponential energy weights
	\eqn{E(h(x)) = } area of the subset of the window to the left of
	the line \eqn{X=x}.
      \item
	if \code{cumulative=FALSE} then 
	\eqn{h(x)} is the marginal integral of 
	the smoothed residual field (see the case \code{which="smooth"}
	described above) on the \eqn{x} axis. 
	This is approximately the derivative
	of the plot for \code{cumulative=TRUE}.
	The value of \eqn{h(x)} is computed by summing the values of
	the smoothed residual field over all pixels with
	the given \eqn{x} coordinate. 
	For the residuals \eqn{E(h(x)) = 0},
	and for the exponential energy weights
	\eqn{E(h(x)) = } length of the intersection between the
	observation window and the line \eqn{X=x}.
      }
      If \code{plot.sd = TRUE}, then superimposed on the lurking variable
      plot are the pointwise
      two-standard-deviation error limits for \eqn{h(x)} calculated for the
      inhomogeneous Poisson process. The default is \code{plot.sd = TRUE}
      for Poisson models and \code{plot.sd = FALSE} for non-Poisson
      models.
    }
    \item{\code{"y"}:}{
      produce a similar lurking variable plot for the \eqn{y} coordinate.
    }
    \item{\code{"sum"}:}{
      print the sum of the weights or residuals for all points
      in the window (clipped by a margin \code{rbord} if required)
      and the area of the same window. If the fitted model is correct
      the sum of the exponential energy weights should equal the area of
      the window, while the sum of the residuals should equal zero.
      Also print the range of values of the smoothed field
      displayed in the \code{"smooth"} case.
    }
    \item{\code{"all"}:}{
      All four of the diagnostic plots listed above are plotted together
      in a two-by-two display. Top left panel is \code{"marks"} plot.
      Bottom right panel is \code{"smooth"} plot. Bottom left panel is
      \code{"x"} plot. Top right panel is \code{"y"} plot, rotated 90 degrees.
    }
  }

  The argument \code{rbord} ensures there are no edge
  effects in the computation of the residuals. The diagnostic calculations
    will be confined to those points of the data pattern which are
    at least \code{rbord} units away from the edge of the window.
  The value of \code{rbord} should be greater than or equal to
  the range of interaction permitted in the model.

  The argument \code{rv} would normally be used only by experts.
  It enables the user to substitute arbitrary values for the
  residuals or marks, overriding the usual calculations.
  If \code{rv} is present, then instead of calculating the residuals from
  the fitted model, the algorithm takes the entries of the vector
  \code{rv} to be the values of the residuals or marks,
  and plots them in the manner appropriate to the type of residual
  or mark selected by \code{type}. The length of \code{rv} must equal
  the number of points in the original data point pattern
  (if \code{type="eem"}) or the number of quadrature (data and dummy)
  points used to fit the model (if \code{type="raw"},
  \code{type="inverse"} or \code{type="pearson"}).

  The return value is an object of class \code{"diagppm"}.
  There are methods for \code{plot} and \code{print} for such objects.
  See the Examples.
  
  See also the companion function \code{\link{qqplot.ppm}} which produces a
  Q-Q plot of the residuals.

}
\references{
  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  Journal of the Royal Statistical Society, Series B.
  to appear.

  Stoyan, D. and Grabarnik, P. (1991)
  Second-order characteristics for stochastic structures connected with
  Gibbs point processes.
  \emph{Mathematische Nachrichten}, 151:95--100.
}
\seealso{
 \code{\link{residuals.ppm}},
 \code{\link{eem}},
 \code{\link{ppm.object}},
 \code{\link{qqplot.ppm}},
 \code{\link{ppm}}
}
\examples{
    data(cells)
    fit <- ppm(cells, ~x, Strauss(r=0.15), rbord=0.15)
    diagnose.ppm(fit, rbord=0.15)
    diagnose.ppm(fit, type="pearson", rbord=0.15)
    diagnose.ppm(fit, rbord=0.15, which="marks")
    diagnose.ppm(fit, rbord=0.15, type="raw", plot.neg="discrete")

    # save the diagnostics and plot them later
    u <- diagnose.ppm(fit, rbord=0.15, plot.it=FALSE)
    plot(u)
    plot(u, which="marks") 
}
\author{Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\keyword{spatial}
