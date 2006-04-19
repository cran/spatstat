\name{spatstat}
\alias{spatstat}
\alias{spatstat-package}
\docType{package}
\title{
The Spatstat Package
}
\description{
  This is a summary of the features of 
  \pkg{spatstat}, a package in \code{R}
  for the statistical analysis of spatial point patterns.
}
\details{
  \pkg{spatstat} is a package for the statistical analysis
  of spatial data. Currently, it deals mainly with the analysis of
  patterns of points in the plane.  The points may carry `marks',
  and the spatial region in which the points were recorded 
  may have arbitrary shape. 

  The package supports
  \itemize{
    \item creation, manipulation and plotting of point patterns
    \item exploratory data analysis 
    \item simulation of point process models
    \item parametric model-fitting
    \item hypothesis tests and diagnostics
  }
  The point process models to be fitted
  may be quite general Gibbs/Markov models; they may include spatial trend,
  dependence on covariates, and interpoint interactions of any order (i.e.
  not restricted to pairwise interactions). Models are specified by 
  a \code{formula} in the \code{R} language, and are fitted using
  a single function \code{\link{ppm}} analogous to
  \code{\link{lm}} and \code{\link{glm}}.

}
\section{Getting Started}{
  Type \code{demo(spatstat)} for an overall demonstration 
  of the package.

  An introduction to the package is available in the
  file \code{Intro.pdf} in the \code{doc} directory supplied in the
  installation of the package. 

  The package is described more fully in the paper by
  Baddeley and Turner (2005a), available online.
}
\section{FUNCTIONS AND DATASETS}{
  Following is a summary of the main functions and datasets
  in the \pkg{spatstat} package.
  Alternatively an alphabetical list of all functions and
  datasets is available by typing \code{library(help=spatstat)}.

  For further information on any of these,
  type \code{help(name)} where \code{name} is the name of the function
  or dataset.
}
\section{I. CREATING AND MANIPULATING DATA}{

  \bold{To create a point pattern:}
  
  \tabular{ll}{
    \code{\link{ppp}} \tab
    create a point pattern from \eqn{(x,y)} and window information
    \cr
    \tab
    \code{ppp(x, y, xlim, ylim)} for rectangular window\cr
    \tab
    \code{ppp(x, y, poly)} for polygonal window \cr
    \tab
    \code{ppp(x, y, mask)} for binary image window \cr
    \code{\link{as.ppp}} \tab
    convert other types of data to a \code{ppp} object \cr
    \code{\link{setmarks}}, \code{\%mark\%}  \tab
    attach/reassign marks to a point pattern
  }
      
  \bold{To simulate a random point pattern:}
  
  \tabular{ll}{
	\code{\link{runifpoint}} \tab
	generate \eqn{n} independent uniform random points \cr
	\code{\link{rpoint}} \tab
	generate \eqn{n} independent random points \cr
	\code{\link{rmpoint}} \tab
	generate \eqn{n} independent multitype random points \cr
	\code{\link{rpoispp}} \tab
	simulate the (in)homogeneous Poisson point process \cr
	\code{\link{rmpoispp}} \tab
	simulate the (in)homogeneous multitype Poisson point process \cr
	\code{\link{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link{rsyst}} \tab
	systematic random sample of points \cr
	\code{\link{rMaternI}}  \tab
	simulate the Mat\'ern Model I inhibition process\cr
	\code{\link{rMaternII}} \tab
	simulate the Mat\'ern Model II inhibition process\cr
	\code{\link{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link{rMatClust}} \tab
	simulate the Mat\'ern Cluster process\cr
	\code{\link{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
	\code{\link{rshift}} \tab random shift (including toroidal shifts)
      }

      \bold{Standard point pattern datasets:}
      
      Remember to say \code{\link{data}(bramblecanes)} etc.
      
      \tabular{ll}{
	\code{\link{amacrine}} \tab Austin Hughes' rabbit amacrine cells \cr
	\code{\link{ants}} \tab Harkness-Isham ant nests data\cr
	\code{\link{betacells}} \tab Waessle et al. cat retinal ganglia data \cr
	\code{\link{bramblecanes}} \tab Bramble Canes data \cr
	\code{\link{cells}} \tab Crick-Ripley biological cells data \cr
	\code{\link{chorley}} \tab Chorley-Ribble cancer data \cr
	\code{\link{copper}} \tab Berman-Huntington copper deposits data \cr
	\code{\link{demopat}} \tab Synthetic point pattern \cr
	\code{\link{finpines}} \tab Finnish Pines data \cr
	\code{\link{hamster}} \tab Aherne's hamster tumour data \cr
	\code{\link{humberside}} \tab North Humberside childhood leukaemia data \cr
	\code{\link{japanesepines}} \tab Japanese Pines data \cr
	\code{\link{lansing}} \tab Lansing Woods data \cr
	\code{\link{longleaf}} \tab Longleaf Pines data \cr
	\code{\link{nztrees}} \tab Mark-Esler-Ripley trees data \cr
	\code{\link{redwood}} \tab Strauss-Ripley redwood saplings data \cr
	\code{\link{redwoodfull}} \tab Strauss redwood saplings data (full set) \cr
	\code{\link{residualspaper}} \tab Data from Baddeley et al (2005) \cr
	\code{\link{simdat}} \tab Simulated point pattern (inhomogeneous, with interaction) \cr
	\code{\link{spruces}} \tab Spruce trees in Saxonia \cr
	\code{\link{swedishpines}} \tab Strand-Ripley swedish pines data
      }

      \bold{To manipulate a point pattern:}

      \tabular{ll}{
	\code{\link{plot.ppp}} \tab
	plot a point pattern (e.g. \code{plot(X)}) \cr
	\code{\link{subset.ppp}},
	\code{"[.ppp"} \tab
	extract or replace a subset of a point pattern \cr
        \tab \code{pp[subset]} \cr
	\tab  \code{pp[, subwindow]} \cr
	\code{\link{superimpose}} \tab
	superimpose any number of point patterns  \cr
	\code{\link{cut.ppp}} \tab
	discretise the marks in a point pattern \cr
	\code{\link{unmark}} \tab
	remove marks  \cr
	\code{\link{setmarks}} \tab
	attach marks or reset marks  \cr
	\code{\link{split.ppp}} \tab
	divide pattern into sub-patterns \cr
	\code{\link{rotate}} \tab
	rotate pattern  \cr
	\code{\link{shift} } \tab
	translate pattern  \cr
	\code{\link{affine}} \tab
	apply affine transformation\cr
	\code{\link{density.ppp}} \tab
	kernel smoothing\cr
	\code{\link{identify.ppp}} \tab
	interactively identify points \cr
	\code{\link{unique.ppp}} \tab
	remove duplicate points \cr
	\code{\link{duplicated.ppp}} \tab
	determine which points are duplicates
      }
      See \code{\link{spatstat.options}} to control plotting behaviour.
      
    \bold{To create a window:}

    An object of class \code{"owin"} describes a spatial region
      (a window of observation).

      \tabular{ll}{
	\code{\link{owin}}	\tab	Create a window object \cr
	\tab \code{owin(xlim, ylim)} for rectangular window \cr
	\tab \code{owin(poly)} for polygonal window \cr
	\tab \code{owin(mask)} for binary image window \cr
	\code{\link{as.owin}}	\tab
	Convert other data to a window object \cr
	\code{\link{ripras}}    \tab
	Ripley-Rasson estimator of window, given only the points \cr
	\code{\link{letterR}}    \tab
	polygonal window in the shape of the {\sf R} logo
      }

    \bold{To manipulate a window:}

    \tabular{ll}{
	\code{\link{plot.owin}}	\tab	plot a window. \cr
	\tab		\code{plot(W)}\cr
	\code{\link{bounding.box}} \tab
	Find a tight bounding box for the window \cr
	\code{\link{erode.owin}}	\tab
	erode window by a distance r\cr
	\code{\link{dilate.owin}}	\tab
	dilate window by a distance r\cr
	\code{\link{complement.owin}}	\tab
	invert (swap inside and outside)\cr
	\code{\link{rotate}} \tab rotate window  \cr
	\code{\link{shift} } \tab translate window  \cr
	\code{\link{affine}} \tab apply affine transformation
      }

    \bold{Digital approximations:}

    \tabular{ll}{
	\code{\link{as.mask}}	\tab
	Make a discrete pixel approximation of a given window \cr
	\code{\link{nearest.raster.point}} \tab
	map continuous coordinates to raster locations\cr
	\code{\link{raster.x}} \tab
	raster x coordinates \cr
	\code{\link{raster.y}} \tab
	raster y coordinates 
      }
      See \code{\link{spatstat.options}} to control the approximation

    \bold{Geometrical computations with windows:}

    \tabular{ll}{
	\code{\link{intersect.owin}}	\tab	intersection of two windows\cr
	\code{\link{union.owin}}	\tab	union of two windows\cr
	\code{\link{inside.owin}}	\tab	determine whether a point is inside a window\cr
	\code{\link{area.owin}}	\tab	compute window's area \cr
	\code{\link{diameter}}	\tab	compute window frame's diameter\cr
	\code{\link{eroded.areas}}	\tab	compute areas of eroded windows\cr
	\code{\link{bdist.points}}	\tab	compute distances from data points to window boundary \cr
	\code{\link{bdist.pixels}}	\tab	compute distances from all pixels to window boundary \cr
	\code{\link{distmap.owin}}	\tab	distance transform image \cr
	\code{\link{centroid.owin}}	\tab compute centroid (centre of mass) of window\cr
	\code{\link{is.subset.owin}}    \tab determine whether one window contains another 
      }

    \bold{Pixel images:}
    An object of class \code{"im"} represents a pixel image. 
    Such objects are returned by some of the functions in
      \pkg{spatstat} including \code{\link{Kmeasure}},
      \code{\link{setcov}} and \code{\link{density.ppp}}. 
      \tabular{ll}{
	\code{\link{im}} \tab create a pixel image\cr
	\code{\link{as.im}} \tab convert other data to a pixel image\cr
	\code{\link{as.matrix.im}} \tab convert pixel image to matrix\cr
	\code{\link{plot.im}}	\tab	plot a pixel image on screen as a digital image\cr
	\code{\link{contour.im}}	\tab draw contours of a pixel image \cr
	\code{\link{persp.im}}	\tab draw perspective plot of a pixel image \cr
	\code{\link{[.im}} 	\tab extract subset of pixel image\cr
	\code{\link{shift.im}} \tab apply vector shift to pixel image \cr
	\code{X}	\tab	print very basic information about image \code{X}\cr
	\code{\link{summary}(X)} \tab	summary of image \code{X} \cr
	\code{\link{hist.im}} \tab	histogram of image \cr
	\code{\link{mean.im}} \tab	mean pixel value of image  \cr
	\code{\link{quantile.im}} \tab	quantiles of image  \cr
	\code{\link{cut.im}} \tab	convert numeric image to factor image \cr
	\code{\link{is.im}} \tab test whether an object is a pixel image\cr
	\code{\link{compatible.im}} \tab test whether two images have
	compatible dimensions \cr
	\code{\link{eval.im}} \tab evaluate any expression involving images\cr
	\code{\link{levelset}} \tab level set of an image\cr
	\code{\link{solutionset}} \tab region where an expression is true
      }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of line
    segments.
    \tabular{ll}{
      \code{\link{psp}} \tab create a line segment pattern \cr
      \code{\link{as.psp}} \tab convert other data into a line segment pattern \cr
      \code{\link{plot.psp}} \tab plot a line segment pattern \cr
      \code{\link{print.psp}} \tab print basic information \cr
      \code{\link{summary.psp}} \tab print summary information \cr
      \code{\link{subset.psp}} \tab \cr
      \code{\link{[.psp}}
                  \tab extract a subset of a line segment pattern \cr
      \code{\link{midpoints.psp}} \tab compute the midpoints of line segments \cr
      \code{\link{endpoints.psp}} \tab extract the endpoints of line segments \cr
      \code{\link{lengths.psp}} \tab compute the lengths of line segments \cr
      \code{\link{angles.psp}} \tab compute the orientation angles of line segments \cr
      \code{\link{rotate.psp}} \tab rotate a line segment pattern \cr      
      \code{\link{shift.psp}} \tab shift a line segment pattern \cr      
      \code{\link{affine.psp}} \tab apply an affine transformation \cr      
      \code{\link{distmap.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{density.psp}} \tab kernel smoothing of line segments
    }

  }
  
\section{II. EXPLORATORY DATA ANALYSIS}{

  \bold{Inspection of data:}
  \tabular{ll}{
    \code{\link{summary}(X)} \tab
    print useful summary of point pattern \code{X}\cr
    \code{X} \tab print basic description of point pattern \code{X} 
  }

  \bold{Summary statistics for a point pattern:}
  \tabular{ll}{
    \code{\link{quadratcount}} \tab Quadrat counts \cr
    \code{\link{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link{pcf}} \tab 	pair correlation function \cr
    \code{\link{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link{envelope}} \tab simulation envelopes for a summary
    function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link{plot.fv}} \tab plot a summary function\cr
    \code{\link{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link{nndist}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich}} \tab find nearest neighbours \cr
    \code{\link{pairdist}} \tab distances between all pairs of points\cr
    \code{\link{crossdist}} \tab distances between points in two patterns\cr
    \code{\link{exactdt}} \tab distance from any location to nearest data point\cr 
    \code{\link{distmap}} \tab distance map image\cr
    \code{\link{density.ppp}} \tab kernel smoothed density
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks} which is a factor. 
  \tabular{ll}{
    \code{\link{Gcross},\link{Gdot},\link{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{Gij, Gi.} \cr
    \code{\link{Kcross},\link{Kdot}, \link{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{Kij, Ki.} \cr
    \code{\link{Jcross},\link{Jdot},\link{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{Jij,Ji.} \cr
    \code{\link{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link{Iest}} \tab  multitype \eqn{I}-function
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type.
  \tabular{ll}{
    \code{\link{markcorr}} \tab mark correlation function \cr
    \code{\link{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link{applynbd}} \tab apply function to every neighbourhood\cr
    \tab in a point pattern 
  }
}
\section{III. MODEL FITTING}{
  \bold{To fit a point process model:}

  Model fitting in \pkg{spatstat} is performed by the function
  \code{\link{ppm}}. Its result is an object of class \code{"ppm"}. 

  \bold{Manipulating the fitted model:}

  \tabular{ll}{
    \code{\link{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link{predict.ppm}} \tab 		Compute the spatial trend \cr
    \tab	and conditional intensity\cr
    \tab        of the fitted point process model \cr
    \code{\link{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link{update.ppm}} \tab Update the fit \cr
    \code{\link{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link{anova.ppm}} \tab Analysis of deviance
  }
  See \code{\link{spatstat.options}} to control plotting of fitted model.

  
  \bold{To specify a point process model:}
  
  The first order ``trend'' of the model is written as an \code{R}
  language formula.
  
  \tabular{ll}{
    \code{~1}  \tab No trend (stationary) \cr
    \code{~x}   \tab First order term  \eqn{\lambda(x,y) = \exp(\alpha + \beta x)} \cr
    \tab where \eqn{x,y} are Cartesian coordinates \cr
    \code{~polynom(x,y,3)}  \tab Log-cubic polynomial trend  \cr
    \code{~harmonic(x,y,2)}  \tab Log-harmonic polynomial trend 
  }

  The higher order (``interaction'') components are described by
  an object of class \code{"interact"}. Such objects are created by:
  \tabular{ll}{
    \code{\link{Poisson}()} \tab		the Poisson point process\cr
    \code{\link{Strauss}()} \tab		the Strauss process \cr
    \code{\link{StraussHard}()} \tab		the Strauss/hard core point process\cr
    \code{\link{Softcore}()} \tab		pairwise interaction, soft core potential\cr
    \code{\link{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link{MultiStraussHard}()} \tab 	multitype Strauss/hard core process 
  }
  
  \bold{Finer control over model fitting:}
  
  A quadrature scheme is represented by an object of
  class \code{"quad"}.

  \tabular{ll}{
    \code{\link{quadscheme}} \tab
    generate a Berman-Turner quadrature scheme\cr
    \tab	for use by \code{ppm } \cr
    \code{\link{default.dummy}} \tab default pattern of dummy points \cr
    \code{\link{gridcentres}} \tab dummy points in a rectangular grid \cr
    \code{\link{rstrat}} \tab stratified random dummy pattern \cr
    \code{\link{spokes}} \tab radial pattern of dummy points  \cr
    \code{\link{corners}} \tab dummy points at corners of the window  \cr
    \code{\link{gridweights}} \tab quadrature weights by the grid-counting rule  \cr
    \code{\link{dirichlet.weights}} \tab quadrature weights are Dirichlet tile areas \cr
    \code{plot(Q)} \tab plot quadrature scheme \code{Q}\cr
    \code{print(Q)} \tab print basic information about quadrature scheme \code{Q}\cr
    \code{\link{summary}(Q)} \tab summary of quadrature scheme \code{Q}
  }

  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }
}
\section{IV. TESTS AND DIAGNOSTICS}{

  \bold{Classical hypothesis tests:}
  
  \tabular{ll}{
    \code{\link{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link{anova.ppm}} \tab Analysis of Deviance for
    point process models 
  }
  
  \bold{Diagnostic plots:}
  
  Residuals for a fitted point process model, and diagnostic plots
  based on the residuals, were introduced in Baddeley et al (2005).
  
  Type \code{demo(diagnose)}
  for a demonstration of the diagnostics features.

  \tabular{ll}{
    \code{\link{diagnose.ppm}} \tab diagnostic plots for spatial trend\cr
    \code{\link{qqplot.ppm}} \tab diagnostic plot for interpoint interaction\cr
    \code{\link{residualspaper}} \tab examples from Baddeley et al (2005)
  }
}


\section{Documentation}{
  The online manual entries are quite detailed and should be consulted
  first for information about a particular function.
  
  The paper by Baddeley and Turner (2005a) describes the package.
  Baddeley and Turner (2005b) is a more detailed explanation of
  how to fit point process models to data.
}
\references{
  Baddeley, A. and Turner, R. (2005a)
  Spatstat: an R package for analyzing spatial point patterns.
  \emph{Journal of Statistical Software} \bold{12}:6, 1--42.
  URL: \code{www.jstatsoft.org}, ISSN: 1548-7660.

  Baddeley, A. and Turner, R. (2005b)
  Modelling spatial point patterns in R.
  In: A. Baddeley, P. Gregori, J. Mateu, R. Stoica, and D. Stoyan,
  editors, \emph{Case Studies in Spatial Point Pattern Modelling},
  Lecture Notes in Statistics number 185. Pages 23--74.
  Springer-Verlag, New York, 2006. 
  ISBN: 0-387-28311-0.

  Baddeley, A., Turner, R., Moller, J. and Hazelton, M. (2005)
  Residual analysis for spatial point processes.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{67}, 617--666.
}
\section{Licence}{
  This library and its documentation are usable under the terms of the "GNU 
  General Public License", a copy of which is distributed with the package.
}
\author{
  Adrian Baddeley
  \email{adrian@maths.uwa.edu.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{rolf@math.unb.ca}
  \url{http://www.math.unb.ca/~rolf}
}
\section{Acknowledgements}{
  Marie-Colette van Lieshout made substantial contributions of code.
  Contributions were also made by J.B. Chen, Y.C. Chin, S. Eglen,
  P. Grabarnik, U. Hahn, M.B. Hansen, K. Hornik, J. Mateu, L.S. Nielsen,
  B.D. Ripley, B. Rowlingson, A. Sarkka, K. Schladitz,
  B.T. Scott, M. Stevenson, B. Turlach, and A. van Burgel.
}
\keyword{spatial}

