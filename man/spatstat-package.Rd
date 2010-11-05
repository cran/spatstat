\name{spatstat-package}
\alias{spatstat-package}
\alias{spatstat}
\docType{package}
\title{The Spatstat Package}
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
  It is also possible to fit cluster process models by the method of
  minimum contrast.

  Apart from two-dimensional point patterns and point processes,
  \pkg{spatstat} also supports patterns of line segments in two
  dimensions, point patterns in three dimensions, and 
  multidimensional space-time point patterns. It also supports
  spatial tessellations and random sets.
}
\section{Getting Started}{
  For a quick introduction to \pkg{spatstat}, see 
  the paper by Baddeley and Turner (2005a).
  For a complete 2-day course on using \code{spatstat}, see the workshop notes
  by Baddeley (2008).
  Both of these documents are available on the internet.

  Type \code{demo(spatstat)} for a demonstration 
  of the package's capabilities.
  Type \code{demo(data)} to see all the datasets
  available in the package.
}
\section{Updates}{
  New versions of \pkg{spatstat} are produced about once a month.
  Users are advised to update their installation of \pkg{spatstat}
  regularly.
  
  Type \code{latest.news()} to read the news documentation about
  changes to the current installed version of \pkg{spatstat}. Type
  \code{news(package="spatstat")} to read news documentation about
  all previous versions of the package.
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
\section{CONTENTS:}{
  \tabular{ll}{
    I. \tab Creating and manipulating data \cr
    II. \tab Exploratory Data Analysis \cr
    III. \tab Model fitting (cluster models) \cr
    IV. \tab Model fitting (Poisson and Gibbs models) \cr
    V. \tab Model fitting (spatial logistic regression)\cr
    VI. \tab Simulation \cr
    VII. \tab Tests and diagnostics\cr
    VIII. \tab Documentation
  }
}

\section{I. CREATING AND MANIPULATING DATA}{

  \bold{Types of spatial data:}

  The main types of spatial data supported by \pkg{spatstat} are:

  \tabular{ll}{
    \code{\link{ppp}} \tab point pattern \cr
    \code{\link{owin}} \tab window (spatial region) \cr
    \code{\link{im}} \tab pixel image \cr
    \code{\link{psp}} \tab line segment pattern \cr
    \code{\link{tess}} \tab tessellation \cr
    \code{\link{pp3}} \tab three-dimensional point pattern \cr
    \code{\link{ppx}} \tab point pattern in any number of dimensions
  }

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
    \code{\link{clickppp}} \tab
    interactively add points to a plot \cr
    \code{\link{marks<-}}, \code{\%mark\%}  \tab
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
	\code{\link{runifdisc}} \tab
	generate \eqn{n} independent uniform random points in disc\cr
	\code{\link{rstrat}} \tab
	stratified random sample of points \cr
	\code{\link{rsyst}} \tab
	systematic random sample of points \cr
	\code{\link{rjitter}} \tab
	apply random displacements to points in a pattern\cr
	\code{\link{rMaternI}}  \tab
	simulate the Mat\'ern Model I inhibition process\cr
	\code{\link{rMaternII}} \tab
	simulate the Mat\'ern Model II inhibition process\cr
	\code{\link{rSSI}} \tab
	simulate Simple Sequential Inhibition process\cr
	\code{\link{rStrauss}} \tab
	simulate Strauss process (perfect simulation)\cr
	\code{\link{rNeymanScott}} \tab
	simulate a general Neyman-Scott process\cr
	\code{\link{rMatClust}} \tab
	simulate the Mat\'ern Cluster process\cr
	\code{\link{rThomas}} \tab
	simulate the Thomas process  \cr
	\code{\link{rGaussPoisson}}  \tab
	simulate the Gauss-Poisson cluster process\cr
	\code{\link{rthin}} \tab  random thinning  \cr
	\code{\link{rcell}} \tab
	simulate the Baddeley-Silverman cell process  \cr
	\code{\link{rmh}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link{simulate.ppm}} \tab
	simulate Gibbs point process using Metropolis-Hastings \cr
	\code{\link{runifpointOnLines}} \tab
	generate \eqn{n} random points along specified line segments \cr
	\code{\link{rpoisppOnLines}} \tab
	generate Poisson random points along specified line segments 
      }

      \bold{To randomly change an existing point pattern:}
      
      \tabular{ll}{
	\code{\link{rlabel}} \tab random (re)labelling of a multitype
	point pattern \cr
	\code{\link{rshift}} \tab random shift (including toroidal shifts)
      }

      \bold{Standard point pattern datasets:}
      
      Remember to say \code{\link{data}(amacrine)} etc.
      Type \code{demo(data)} to see all the datasets
      installed with the package.
      
      \tabular{ll}{
	\code{\link{amacrine}} \tab Austin Hughes' rabbit amacrine cells \cr
	\code{\link{anemones}} \tab Upton-Fingleton sea anemones data\cr
	\code{\link{ants}} \tab Harkness-Isham ant nests data\cr
	\code{\link{bei}} \tab Tropical rainforest trees \cr
	\code{\link{betacells}} \tab Waessle et al. cat retinal ganglia data \cr
	\code{\link{bramblecanes}} \tab Bramble Canes data \cr
	\code{\link{bronzefilter}} \tab Bronze Filter Section data \cr
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
	\code{\link{murchison}} \tab Murchison gold deposits \cr
	\code{\link{nbfires}} \tab New Brunswick fires data \cr
	\code{\link{nztrees}} \tab Mark-Esler-Ripley trees data \cr
	\code{\link{osteo}} \tab Osteocyte lacunae (3D, replicated) \cr
	\code{\link{ponderosa}} \tab Getis-Franklin ponderosa pine trees data \cr
	\code{\link{redwood}} \tab Strauss-Ripley redwood saplings data \cr
	\code{\link{redwoodfull}} \tab Strauss redwood saplings data (full set) \cr
	\code{\link{residualspaper}} \tab Data from Baddeley et al (2005) \cr
	\code{\link{shapley}} \tab Galaxies in an astronomical survey \cr
	\code{\link{simdat}} \tab Simulated point pattern (inhomogeneous, with interaction) \cr
	\code{\link{spruces}} \tab Spruce trees in Saxonia \cr
	\code{\link{swedishpines}} \tab Strand-Ripley swedish pines data \cr
	\code{\link{urkiola}} \tab Urkiola Woods data
      }

      \bold{To manipulate a point pattern:}

      \tabular{ll}{
	\code{\link{plot.ppp}} \tab
	plot a point pattern (e.g. \code{plot(X)}) \cr
	\code{\link{iplot}} \tab
	plot a point pattern interactively \cr
	\code{\link{[.ppp}} \tab
	extract or replace a subset of a point pattern \cr
        \tab \code{pp[subset]} or \code{pp[subwindow]} \cr
	\code{\link{superimpose}} \tab
	combine several point patterns  \cr
	\code{\link{by.ppp}} \tab
	apply a function to sub-patterns of a point pattern \cr
	\code{\link{cut.ppp}} \tab
	classify the points in a point pattern \cr
	\code{\link{unmark}} \tab
	remove marks  \cr
	\code{\link{marks}} \tab
	attach marks or reset marks  \cr
	\code{\link{npoints}} \tab
	count the number of points  \cr
	\code{\link{split.ppp}} \tab
	divide pattern into sub-patterns \cr
	\code{\link{rotate}} \tab
	rotate pattern  \cr
	\code{\link{shift} } \tab
	translate pattern  \cr
	\code{\link{periodify} } \tab
	make several translated copies  \cr
	\code{\link{affine}} \tab
	apply affine transformation\cr
	\code{\link{density.ppp}} \tab
	kernel smoothing\cr
	\code{\link{smooth.ppp}} \tab
	interpolate pixel values\cr
	\code{\link{sharpen.ppp}} \tab
	data sharpening\cr
	\code{\link{identify.ppp}} \tab
	interactively identify points \cr
	\code{\link{unique.ppp}} \tab
	remove duplicate points \cr
	\code{\link{duplicated.ppp}} \tab
	determine which points are duplicates \cr
	\code{\link{dirichlet}} \tab
	compute Dirichlet-Voronoi tessellation \cr
	\code{\link{delaunay}} \tab
	compute Delaunay triangulation \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{discretise}} \tab discretise coordinates \cr
	\code{\link{pixellate.ppp}} \tab approximate point pattern by 
	pixel image \cr
	\code{\link{as.im.ppp}} \tab approximate point pattern by 
	pixel image 
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
	\code{\link{square}}    \tab	make a square window \cr
	\code{\link{disc}}    \tab	make a circular window \cr
	\code{\link{ripras}}    \tab
	Ripley-Rasson estimator of window, given only the points \cr
	\code{\link{convexhull}} \tab compute convex hull of something \cr 
	\code{\link{letterR}}    \tab
	polygonal window in the shape of the \R logo
      }

    \bold{To manipulate a window:}

    \tabular{ll}{
	\code{\link{plot.owin}}	\tab	plot a window. \cr
	\tab		\code{plot(W)}\cr
	\code{\link{bounding.box}} \tab
	Find a tight bounding box for the window \cr
	\code{\link{erosion}}	\tab
	erode window by a distance r\cr
	\code{\link{dilation}}	\tab
	dilate window by a distance r\cr
	\code{\link{closing}}	\tab
	close window by a distance r\cr
	\code{\link{opening}}	\tab
	open window by a distance r\cr
	\code{\link{border}}	\tab
	difference between window and its erosion/dilation \cr 
	\code{\link{complement.owin}}	\tab
	invert (swap inside and outside)\cr
	\code{\link{simplify.owin}}	\tab
	approximate a window by a simple polygon  \cr
	\code{\link{rotate}} \tab rotate window  \cr
	\code{\link{shift} } \tab translate window  \cr
	\code{\link{periodify} } \tab make several translated copies  \cr
	\code{\link{affine}} \tab apply affine transformation 
      }

    \bold{Digital approximations:}

    \tabular{ll}{
	\code{\link{as.mask}}	\tab
	Make a discrete pixel approximation of a given window \cr
	\code{\link{as.im.owin}} \tab convert window to pixel image \cr
	\code{\link{pixellate.owin}} \tab convert window to pixel image \cr
	\code{\link{nearest.raster.point}} \tab
	map continuous coordinates to raster locations\cr
	\code{\link{raster.x}} \tab
	raster x coordinates \cr
	\code{\link{raster.y}} \tab
	raster y coordinates \cr
	\code{\link{as.polygonal}} \tab
	convert pixel mask to polygonal window
      }
      See \code{\link{spatstat.options}} to control the approximation

    \bold{Geometrical computations with windows:}

    \tabular{ll}{
	\code{\link{intersect.owin}}	\tab	intersection of two windows\cr
	\code{\link{union.owin}}	\tab	union of two windows\cr
	\code{\link{setminus.owin}}	\tab	set subtraction of two windows\cr
	\code{\link{inside.owin}}	\tab	determine whether a point is inside a window\cr
	\code{\link{area.owin}}	\tab	compute area \cr
	\code{\link{perimeter}}	\tab	compute perimeter length \cr
	\code{\link{diameter.owin}}	\tab	compute diameter\cr
	\code{\link{incircle}}	\tab	find largest circle inside a window \cr
	\code{\link{connected}}    \tab find connected components of window \cr
	\code{\link{eroded.areas}}	\tab	compute areas of eroded windows\cr
	\code{\link{dilated.areas}}	\tab	compute areas of dilated windows\cr
	\code{\link{bdist.points}}	\tab	compute distances from data points to window boundary \cr
	\code{\link{bdist.pixels}}	\tab	compute distances from all pixels to window boundary \cr
	\code{\link{bdist.tiles}}	\tab
	boundary distance for each tile in tessellation \cr
	\code{\link{distmap.owin}}	\tab	distance transform image \cr
	\code{\link{distfun.owin}}	\tab	distance transform \cr
	\code{\link{centroid.owin}}	\tab compute centroid (centre of mass) of window\cr
	\code{\link{is.subset.owin}}    \tab determine whether one
	window contains another \cr
	\code{\link{is.convex}} \tab determine whether a window is convex \cr 
	\code{\link{convexhull}} \tab compute convex hull \cr
	\code{\link{as.mask}} \tab pixel approximation of window \cr
	\code{\link{as.polygonal}} \tab polygonal approximation of window 
      }

    \bold{Pixel images:}
    An object of class \code{"im"} represents a pixel image. 
    Such objects are returned by some of the functions in
      \pkg{spatstat} including \code{\link{Kmeasure}},
      \code{\link{setcov}} and \code{\link{density.ppp}}. 
      \tabular{ll}{
	\code{\link{im}} \tab create a pixel image\cr
	\code{\link{as.im}} \tab convert other data to a pixel image\cr
	\code{\link{pixellate}} \tab convert other data to a pixel image\cr
	\code{\link{as.matrix.im}} \tab convert pixel image to matrix\cr
	\code{\link{as.data.frame.im}} \tab convert pixel image to data frame\cr
	\code{\link{plot.im}}	\tab	plot a pixel image on screen as a digital image\cr
	\code{\link{contour.im}}	\tab draw contours of a pixel image \cr
	\code{\link{persp.im}}	\tab draw perspective plot of a pixel image \cr
	\code{\link{[.im}} 	\tab extract a subset of a pixel image\cr
	\code{\link{[<-.im}} 	\tab replace a subset of a pixel image\cr
	\code{\link{shift.im}} \tab apply vector shift to pixel image \cr
	\code{X}	\tab	print very basic information about image \code{X}\cr
	\code{\link{summary}(X)} \tab	summary of image \code{X} \cr
	\code{\link{hist.im}} \tab	histogram of image \cr
	\code{\link{mean.im}} \tab	mean pixel value of image  \cr
	\code{\link{integral.im}} \tab	integral of pixel values  \cr
	\code{\link{quantile.im}} \tab	quantiles of image  \cr
	\code{\link{cut.im}} \tab	convert numeric image to factor image \cr
	\code{\link{is.im}} \tab test whether an object is a pixel image\cr
	\code{\link{interp.im}} \tab interpolate a pixel image\cr
	\code{\link{blur}} \tab apply Gaussian blur to image\cr
	\code{\link{connected}} \tab find connected components \cr
	\code{\link{compatible.im}} \tab test whether two images have
	compatible dimensions \cr
	\code{\link{eval.im}} \tab evaluate any expression involving images\cr
	\code{\link{levelset}} \tab level set of an image\cr
	\code{\link{solutionset}} \tab region where an expression is true
      }

    \bold{Line segment patterns}

    An object of class \code{"psp"} represents a pattern of straight line
    segments.
    \tabular{ll}{
      \code{\link{psp}} \tab create a line segment pattern \cr
      \code{\link{as.psp}} \tab convert other data into a line segment pattern \cr
      \code{\link{is.psp}} \tab determine whether a dataset has class \code{"psp"} \cr
      \code{\link{plot.psp}} \tab plot a line segment pattern \cr
      \code{\link{print.psp}} \tab print basic information \cr
      \code{\link{summary.psp}} \tab print summary information \cr
      \code{\link{[.psp}} \tab extract a subset of a line segment pattern \cr
      \code{\link{as.data.frame.psp}} \tab
      convert line segment pattern to data frame \cr
      \code{\link{marks.psp}} \tab extract marks of line segments \cr
      \code{\link{marks<-.psp}} \tab assign new marks to line segments \cr
      \code{\link{unmark.psp}} \tab delete marks from line segments \cr
      \code{\link{midpoints.psp}} \tab
      compute the midpoints of line segments \cr
      \code{\link{endpoints.psp}} \tab extract the endpoints of line segments \cr
      \code{\link{lengths.psp}} \tab compute the lengths of line segments \cr
      \code{\link{angles.psp}} \tab compute the orientation angles of line segments \cr
      \code{\link{rotate.psp}} \tab rotate a line segment pattern \cr      
      \code{\link{shift.psp}} \tab shift a line segment pattern \cr      
      \code{\link{periodify}} \tab make several shifted copies \cr      
      \code{\link{affine.psp}} \tab apply an affine transformation \cr      
      \code{\link{pixellate.psp}} \tab approximate line segment pattern
    by pixel image \cr      
      \code{\link{as.mask.psp}} \tab approximate line segment pattern
    by binary mask \cr      
      \code{\link{distmap.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{distfun.psp}} \tab compute the distance map of a line
      segment pattern \cr
      \code{\link{density.psp}} \tab kernel smoothing of line segments\cr
      \code{\link{selfcrossing.psp}} \tab find crossing points between
      line segments \cr
      \code{\link{crossing.psp}} \tab find crossing points between
      two line segment patterns \cr
      \code{\link{nncross}} \tab find distance to nearest line segment
      from a given point\cr
      \code{\link{nearestsegment}} \tab find line segment closest to a
      given point \cr
      \code{\link{project2segment}} \tab find location along a line segment
      closest to a given point \cr
      \code{\link{pointsOnLines}} \tab generate points evenly spaced
      along line segment \cr
      \code{\link{rpoisline}} \tab generate a realisation of the
      Poisson line process inside a window\cr
      \code{\link{rlinegrid}} \tab generate a random array of parallel
      lines through a window
    }

    \bold{Tessellations}

    An object of class \code{"tess"} represents a tessellation.

    \tabular{ll}{
      \code{\link{tess}} \tab create a tessellation \cr
      \code{\link{quadrats}} \tab create a tessellation of rectangles\cr
      \code{\link{as.tess}} \tab convert other data to a tessellation \cr
      \code{\link{plot.tess}} \tab plot a tessellation \cr
      \code{\link{tiles}} \tab extract all the tiles of a tessellation \cr
      \code{\link{[.tess}} \tab extract some tiles of a tessellation \cr
      \code{\link{[<-.tess}} \tab change some tiles of a tessellation \cr
      \code{\link{intersect.tess}} \tab intersect two tessellations \cr
                              \tab or restrict a tessellation to a window \cr
      \code{\link{chop.tess}} \tab subdivide a tessellation by a line \cr
      \code{\link{dirichlet}} \tab compute Dirichlet-Voronoi tessellation of points\cr
      \code{\link{delaunay}} \tab compute Delaunay triangulation of points\cr
      \code{\link{rpoislinetess}} \tab generate tessellation using Poisson line
      process 
    }


    \bold{Three-dimensional point patterns}

    An object of class \code{"pp3"} represents a three-dimensional
    point pattern in a rectangular box. The box is represented by
    an object of class \code{"box3"}.

    \tabular{ll}{
      \code{\link{pp3}} \tab create a 3-D point pattern \cr
      \code{\link{plot.pp3}} \tab plot a 3-D point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{unitname.pp3}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{runifpoint3}} \tab generate uniform random points in 3-D \cr
      \code{\link{rpoispp3}} \tab generate Poisson random points in 3-D \cr
      \code{\link{envelope.pp3}} \tab generate simulation envelopes for
      3-D pattern \cr
      \code{\link{box3}} \tab create a 3-D rectangular box \cr
      \code{\link{as.box3}} \tab convert data to 3-D rectangular box \cr
      \code{\link{unitname.box3}} \tab name of unit of length \cr
      \code{\link{diameter.box3}} \tab diameter of box \cr
      \code{\link{volume.box3}} \tab volume of box \cr
      \code{\link{shortside.box3}} \tab shortest side of box \cr
      \code{\link{eroded.volumes}} \tab volumes of erosions of box 
    }

    \bold{Multi-dimensional space-time point patterns}

    An object of class \code{"ppx"} represents a 
    point pattern in multi-dimensional space and/or time.

    \tabular{ll}{
      \code{\link{ppx}} \tab create a multidimensional space-time point pattern \cr
      \code{\link{coords}} \tab extract coordinates \cr
      \code{\link{as.hyperframe}} \tab extract coordinates \cr
      \code{\link{unitname.ppx}} \tab name of unit of length \cr
      \code{\link{npoints}} \tab count the number of points  \cr
      \code{\link{runifpointx}} \tab generate uniform random points \cr
      \code{\link{rpoisppx}} \tab generate Poisson random points \cr
      \code{\link{boxx}} \tab define multidimensional box  \cr
      \code{\link{diameter.boxx}} \tab diameter of box \cr
      \code{\link{volume.boxx}} \tab volume of box \cr
      \code{\link{shortside.boxx}} \tab shortest side of box \cr
      \code{\link{eroded.volumes.boxx}} \tab volumes of erosions of box 
    }

    \bold{Hyperframes}

    A hyperframe is like a data frame, except that the entries
    may be objects of any kind.

    \tabular{ll}{
      \code{\link{hyperframe}} \tab create a hyperframe \cr
      \code{\link{as.hyperframe}} \tab convert data to hyperframe \cr
      \code{\link{plot.hyperframe}} \tab plot hyperframe \cr
      \code{\link{with.hyperframe}} \tab evaluate expression using each row
      of hyperframe \cr
      \code{\link{cbind.hyperframe}} \tab combine hyperframes by columns\cr
      \code{\link{rbind.hyperframe}} \tab combine hyperframes by rows\cr
      \code{\link{as.data.frame.hyperframe}} \tab convert hyperframe to
      data frame 
    }
}
  
\section{II. EXPLORATORY DATA ANALYSIS}{

  \bold{Inspection of data:}
  \tabular{ll}{
    \code{\link{summary}(X)} \tab
    print useful summary of point pattern \code{X}\cr
    \code{X} \tab
    print basic description of point pattern \code{X}  \cr
    \code{any(duplicated(X))} \tab
    check for duplicated points in pattern \code{X} \cr
    \code{\link{istat}(X)} \tab Interactive exploratory analysis 
  }

  \bold{Classical exploratory tools:}
  \tabular{ll}{
    \code{\link{clarkevans}} \tab Clark and Evans aggregation index \cr
    \code{\link{fryplot}} \tab Fry plot \cr
    \code{\link{miplot}} \tab Morishita Index plot
  }

  \bold{Modern exploratory tools:}
  \tabular{ll}{
    \code{\link{nnclean}} \tab Byers-Raftery feature detection 
  }

  \bold{Summary statistics for a point pattern:}
  \tabular{ll}{
    \code{\link{quadratcount}} \tab Quadrat counts \cr
    \code{\link{Fest}} \tab empty space function \eqn{F} \cr
    \code{\link{Gest}} \tab nearest neighbour distribution function \eqn{G} \cr
    \code{\link{Kest}} \tab Ripley's \eqn{K}-function\cr
    \code{\link{Lest}} \tab Ripley's \eqn{L}-function\cr
    \code{\link{Jest}} \tab \eqn{J}-function \eqn{J = (1-G)/(1-F)} \cr
    \code{\link{allstats}} \tab all four functions \eqn{F}, \eqn{G}, \eqn{J}, \eqn{K} \cr
    \code{\link{localL}} \tab Getis-Franklin neighbourhood density function\cr
    \code{\link{localK}} \tab neighbourhood K-function\cr
    \code{\link{pcf}} \tab 	pair correlation function \cr
    \code{\link{Kinhom}} \tab \eqn{K} for inhomogeneous point patterns \cr
    \code{\link{Linhom}} \tab \eqn{L} for inhomogeneous point patterns \cr
    \code{\link{pcfinhom}} \tab pair correlation for inhomogeneous patterns\cr
    \code{\link{Kest.fft}} \tab fast \eqn{K}-function using FFT for large datasets \cr
    \code{\link{Kmeasure}} \tab reduced second moment measure \cr
    \code{\link{envelope}} \tab simulation envelopes for a summary
    function \cr
  }

  Related facilities:
  \tabular{ll}{
    \code{\link{plot.fv}} \tab plot a summary function\cr
    \code{\link{eval.fv}} \tab evaluate any expression involving
    summary functions\cr
    \code{\link{eval.fasp}} \tab evaluate any expression involving
    an array of functions\cr
    \code{\link{with.fv}} \tab evaluate an expression for a 
    summary function\cr
    \code{\link{smooth.fv}} \tab apply smoothing to a summary function\cr
    \code{\link{nndist}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich}} \tab find nearest neighbours \cr
    \code{\link{pairdist}} \tab distances between all pairs of points\cr
    \code{\link{crossdist}} \tab distances between points in two patterns\cr
    \code{\link{nncross}} \tab nearest neighbours between two point patterns \cr
    \code{\link{exactdt}} \tab distance from any location to nearest data point\cr 
    \code{\link{distmap}} \tab distance map image\cr
    \code{\link{distfun}} \tab distance map function\cr
    \code{\link{density.ppp}} \tab kernel smoothed density\cr
    \code{\link{smooth.ppp}} \tab spatial interpolation of marks  \cr
    \code{\link{sharpen.ppp}} \tab data sharpening  \cr
    \code{\link{rknn}} \tab theoretical distribution of nearest
    neighbour distance
 }

  \bold{Summary statistics for a multitype point pattern:}
  A multitype point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks} which is a factor. 
  \tabular{ll}{
    \code{\link{Gcross},\link{Gdot},\link{Gmulti}} \tab
    multitype nearest neighbour distributions 
    \eqn{G_{ij}, G_{i\bullet}}{G[i,j], G[i.]} \cr
    \code{\link{Kcross},\link{Kdot}, \link{Kmulti}} \tab
    multitype \eqn{K}-functions 
    \eqn{K_{ij}, K_{i\bullet}}{K[i,j], K[i.]} \cr
    \code{\link{Lcross},\link{Ldot}} \tab
    multitype \eqn{L}-functions 
    \eqn{L_{ij}, L_{i\bullet}}{L[i,j], L[i.]} \cr
    \code{\link{Jcross},\link{Jdot},\link{Jmulti}} \tab
    multitype \eqn{J}-functions
    \eqn{J_{ij}, J_{i\bullet}}{J[i,j],J[i.]} \cr
    \code{\link{pcfcross}} \tab
    multitype pair correlation function \eqn{g_{ij}}{g[i,j]} \cr
    \code{\link{pcfdot}} \tab
    multitype pair correlation function \eqn{g_{i\bullet}}{g[i.]} \cr
    \code{\link{markconnect}} \tab
    marked connection function \eqn{p_{ij}}{p[i,j]} \cr
    \code{\link{alltypes}} \tab  estimates of the above
    for all \eqn{i,j} pairs \cr
    \code{\link{Iest}} \tab  multitype \eqn{I}-function\cr
    \code{\link{Kcross.inhom},\link{Kdot.inhom}} \tab
    inhomogeneous counterparts of \code{Kcross}, \code{Kdot} \cr
    \code{\link{Lcross.inhom},\link{Ldot.inhom}} \tab
    inhomogeneous counterparts of \code{Lcross}, \code{Ldot} \cr
    \code{\link{pcfcross.inhom},\link{pcfdot.inhom}} \tab
    inhomogeneous counterparts of \code{pcfcross}, \code{pcfdot} 
  }

  \bold{Summary statistics for a marked point pattern:}
  A marked point pattern is represented by an object \code{X}
  of class \code{"ppp"} with a component \code{X$marks}.
  The entries in the vector \code{X$marks} may be numeric, complex,
  string or any other atomic type. For numeric marks, there are the
  following functions:
  \tabular{ll}{
    \code{\link{markmean}} \tab smoothed local average of marks \cr
    \code{\link{markvar}} \tab smoothed local variance of marks \cr
    \code{\link{markcorr}} \tab mark correlation function \cr
    \code{\link{markvario}} \tab mark variogram \cr
    \code{\link{markcorrint}} \tab mark correlation integral \cr
    \code{\link{Emark}} \tab mark independence diagnostic \eqn{E(r)} \cr
    \code{\link{Vmark}} \tab mark independence diagnostic \eqn{V(r)} \cr
    \code{\link{nnmean}} \tab nearest neighbour mean index \cr
    \code{\link{nnvario}} \tab nearest neighbour mark variance index 
  }
  For marks of any type, there are the following:
  \tabular{ll}{
    \code{\link{Gmulti}} \tab multitype nearest neighbour distribution \cr
    \code{\link{Kmulti}} \tab multitype \eqn{K}-function \cr
    \code{\link{Jmulti}} \tab multitype \eqn{J}-function 
  }
  Alternatively use \code{\link{cut.ppp}} to convert a marked point pattern
  to a multitype point pattern.

  \bold{Programming tools:}
  \tabular{ll}{
    \code{\link{applynbd}} \tab apply function to every neighbourhood
    in a point pattern \cr
    \code{\link{markstat}} \tab apply function to the marks of neighbours
    in a point pattern \cr
    \code{\link{marktable}} \tab tabulate the marks of neighbours
    in a point pattern \cr
    \code{\link{pppdist}} \tab find the optimal match between two point
    patterns
  }

  \bold{Summary statistics for a three-dimensional point pattern:}

  These are for 3-dimensional point pattern objects (class \code{pp3}).

  \tabular{ll}{
    \code{\link{F3est}} \tab empty space function \eqn{F} \cr
    \code{\link{G3est}} \tab nearest neighbour function \eqn{G} \cr
    \code{\link{K3est}} \tab \eqn{K}-function \cr
    \code{\link{pcf3est}} \tab pair correlation function
  }

  Related facilities:
  \tabular{ll}{
    \code{\link{envelope.pp3}} \tab simulation envelopes \cr
    \code{\link{pairdist.pp3}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.pp3}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.pp3}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.pp3}} \tab find nearest neighbours
  }

  \bold{Computations for multi-dimensional point pattern:}

  These are for multi-dimensional space-time
  point pattern objects (class \code{ppx}).

  \tabular{ll}{
    \code{\link{pairdist.ppx}} \tab distances between all pairs of
    points \cr
    \code{\link{crossdist.ppx}} \tab distances between points in
    two patterns \cr
    \code{\link{nndist.ppx}} \tab nearest neighbour distances \cr
    \code{\link{nnwhich.ppx}} \tab find nearest neighbours
  }

  \bold{Summary statistics for random sets:}
  
  These work for point patterns (class \code{ppp}),
  line segment patterns (class \code{psp})
  or windows (class \code{owin}).
  
  \tabular{ll}{
    \code{\link{Hest}} \tab spherical contact distribution \eqn{H} \cr
    \code{\link{Gfox}} \tab Foxall \eqn{G}-function \cr
    \code{\link{Jfox}} \tab Foxall \eqn{J}-function
  }
  
}

\section{III. MODEL FITTING (CLUSTER MODELS)}{
  
  Cluster process models (with homogeneous or inhomogeneous intensity) 
  can be fitted by the function \code{\link{kppm}}.
  Its result is an object of class \code{"kppm"}.
  The fitted model can be printed, plotted, predicted, simulated
  and updated.

  \tabular{ll}{
    \code{\link{plot.kppm}} \tab  Plot the fitted model\cr
    \code{\link{predict.kppm}} \tab Compute fitted intensity \cr
    \code{\link{update.kppm}} \tab Update the model \cr
    \code{\link{simulate.kppm}} \tab Generate simulated realisations
  }
  
 Lower-level fitting functions include:

  \tabular{ll}{
    \code{\link{thomas.estK}} \tab fit the Thomas process model \cr
    \code{\link{thomas.estpcf}} \tab fit the Thomas process model \cr
    \code{\link{matclust.estK}} \tab fit the Matern Cluster process model \cr
    \code{\link{matclust.estpcf}} \tab fit the Matern Cluster process model \cr
    \code{\link{lgcp.estK}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link{lgcp.estpcf}} \tab fit a log-Gaussian Cox process model\cr
    \code{\link{mincontrast}} \tab low-level algorithm for fitting models
    \cr \tab by the method of minimum contrast 
  }

  The Thomas and Matern models can also be simulated,
  using \code{\link{rThomas}} and \code{\link{rMatClust}} respectively.
}

\section{IV. MODEL FITTING (POISSON AND GIBBS MODELS)}{
  
  \bold{Types of models}
  
  Poisson point processes are the simplest models for point patterns.
  A Poisson model assumes that the points are stochastically
  independent. It may allow the points to have a non-uniform spatial
  density. The special case of a Poisson process with a uniform
  spatial density is often called Complete Spatial Randomness.
  
  Poisson point processes are included in the more general class of Gibbs point
  process models. In a Gibbs model, there is \emph{interaction}
  or dependence between points. Many different types of interaction
  can be specified.
  
  For a detailed explanation of how to fit Poisson or
  Gibbs point process models to point pattern data using \pkg{spatstat},
  see Baddeley and Turner (2005b) or Baddeley (2008).
  
  \bold{To fit a Poison or Gibbs point process model:}

  Model fitting in \pkg{spatstat} is performed mainly by the function
  \code{\link{ppm}}. Its result is an object of class \code{"ppm"}.
  
  Here are some examples, where \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{ppm(X)} \tab Complete Spatial Randomness \cr
    \code{ppm(X, ~1)} \tab Complete Spatial Randomness \cr
    \code{ppm(X, ~x)} \tab Poisson process with \cr
                             \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{ppm(X, ~1, Strauss(0.1))} \tab Stationary Strauss process \cr
    \code{ppm(X, ~x, Strauss(0.1))} \tab Strauss process with \cr
                              \tab conditional intensity loglinear in \eqn{x}
  }
  It is also possible to fit models that depend on
  other covariates.

  \bold{Manipulating the fitted model:}

  \tabular{ll}{
    \code{\link{plot.ppm}} \tab 		Plot the fitted model\cr
    \code{\link{predict.ppm}}
    \tab   Compute the spatial trend and conditional intensity\cr
    \tab   of the fitted point process model \cr
    \code{\link{coef.ppm}} \tab Extract the fitted model coefficients\cr
    \code{\link{formula.ppm}} \tab Extract the trend formula\cr
    \code{\link{fitted.ppm}} \tab Compute fitted conditional intensity at quadrature points \cr
    \code{\link{residuals.ppm}} \tab Compute point process residuals at quadrature points \cr
    \code{\link{update.ppm}} \tab Update the fit \cr
    \code{\link{vcov.ppm}} \tab Variance-covariance matrix of estimates\cr
    \code{\link{rmh.ppm}} \tab Simulate from fitted model  \cr
    \code{\link{simulate.ppm}} \tab Simulate from fitted model  \cr
    \code{\link{print.ppm}} \tab Print basic information about a fitted model\cr
    \code{\link{summary.ppm}} \tab Summarise a fitted model\cr
    \code{\link{effectfun}} \tab Compute the fitted effect of one covariate\cr
    \code{\link{logLik.ppm}} \tab log-likelihood or log-pseudolikelihood\cr
    \code{\link{anova.ppm}} \tab Analysis of deviance
  }
  For model selection, you can also use 
  the generic functions \code{\link{step}}, \code{\link{drop1}} 
  and \code{\link{AIC}} on fitted point process models.
  
  See \code{\link{spatstat.options}} to control plotting of fitted model.
  
  \bold{To specify a point process model:}
  
  The first order ``trend'' of the model is determined by an \code{R}
  language formula. The formula specifies the form of the
  \emph{logarithm} of the trend.
  
  \tabular{ll}{
    \code{~1}  \tab No trend (stationary) \cr
    \code{~x}   \tab Loglinear trend
      \eqn{\lambda(x,y) = \exp(\alpha + \beta x)}{lambda(x,y) =	exp(alpha + beta * x)} \cr
    \tab where \eqn{x,y} are Cartesian coordinates \cr
    \code{~polynom(x,y,3)}  \tab Log-cubic polynomial trend  \cr
    \code{~harmonic(x,y,2)}  \tab Log-harmonic polynomial trend 
  }

  The higher order (``interaction'') components are described by
  an object of class \code{"interact"}. Such objects are created by:
  \tabular{ll}{
    \code{\link{Poisson}()} \tab		the Poisson point process\cr
    \code{\link{AreaInter}()}	 \tab Area-interaction process\cr
    \code{\link{BadGey}()} \tab	multiscale Geyer process\cr
    \code{\link{DiggleGratton}() } \tab Diggle-Gratton potential \cr
    \code{\link{DiggleGatesStibbard}() } \tab Diggle-Gates-Stibbard potential \cr
    \code{\link{Fiksel}()}	 \tab Fiksel pairwise interaction process\cr
    \code{\link{Geyer}()}	 \tab Geyer's saturation process\cr
    \code{\link{Hardcore}()}	 \tab Hard core process\cr
    \code{\link{LennardJones}() } \tab Lennard-Jones potential \cr
    \code{\link{MultiStrauss}()} \tab 		multitype Strauss process \cr
    \code{\link{MultiStraussHard}()} \tab 	multitype Strauss/hard core process \cr
    \code{\link{OrdThresh}()} \tab		Ord process, threshold potential\cr
    \code{\link{Ord}()} \tab 		        Ord model, user-supplied potential \cr
    \code{\link{PairPiece}()} \tab		pairwise interaction, piecewise constant \cr
    \code{\link{Pairwise}()} \tab	pairwise interaction, user-supplied potential\cr
    \code{\link{SatPiece}()} \tab	Saturated pair model, piecewise  constant potential\cr
    \code{\link{Saturated}()} \tab	Saturated pair model, user-supplied potential\cr
    \code{\link{Softcore}()} \tab pairwise interaction, soft core potential\cr
    \code{\link{Strauss}()} \tab  Strauss process \cr
    \code{\link{StraussHard}()} \tab Strauss/hard core point process
  }
  
  \bold{Finer control over model fitting:}
  
  A quadrature scheme is represented by an object of
  class \code{"quad"}. To create a quadrature scheme, typically
  use \code{\link{quadscheme}}.
  
  \tabular{ll}{
    \code{\link{quadscheme}} \tab default quadrature scheme \cr
                             \tab using rectangular cells or Dirichlet cells\cr
    \code{\link{pixelquad}}  \tab quadrature scheme based on image pixels \cr
    \code{\link{quad}}       \tab create an object of class \code{"quad"}
  }
  
  To inspect a quadrature scheme:
  \tabular{ll}{
    \code{plot(Q)} \tab plot quadrature scheme \code{Q}\cr
    \code{print(Q)} \tab print basic information about quadrature scheme \code{Q}\cr
    \code{\link{summary}(Q)} \tab summary of quadrature scheme \code{Q}
  }

  A quadrature scheme consists of data points, dummy points, and
  weights. To generate dummy points:
  \tabular{ll}{
    \code{\link{default.dummy}} \tab default pattern of dummy points \cr
    \code{\link{gridcentres}} \tab dummy points in a rectangular grid \cr
    \code{\link{rstrat}} \tab stratified random dummy pattern \cr
    \code{\link{spokes}} \tab radial pattern of dummy points  \cr
    \code{\link{corners}} \tab dummy points at corners of the window
  }
  
  To compute weights:
  \tabular{ll}{
    \code{\link{gridweights}} \tab quadrature weights by the grid-counting rule  \cr
    \code{\link{dirichlet.weights}} \tab quadrature weights are
    Dirichlet tile areas
  }

  \bold{Simulation and goodness-of-fit for fitted models:}
  
  \tabular{ll}{
    \code{\link{rmh.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link{simulate.ppm}} \tab simulate realisations of a fitted model \cr
    \code{\link{envelope}} \tab compute simulation envelopes for a
    fitted model 
  }
}

\section{V. MODEL FITTING (SPATIAL LOGISTIC REGRESSION)}{
  
  \bold{Logistic regression}
  
  Pixel-based spatial logistic regression is an alternative
  technique for analysing spatial point patterns
  that is widely used in Geographical Information Systems.
  It is approximately equivalent to fitting a Poisson point process
  model.
  
  In pixel-based logistic regression, the spatial domain is
  divided into small pixels, the presence or absence of a
  data point in each pixel is recorded, and logistic regression
  is used to model the presence/absence indicators as a function
  of any covariates.
  
  Facilities for performing spatial logistic regression are
  provided in \pkg{spatstat} for comparison purposes.
  
  \bold{Fitting a spatial logistic regression}
  
  Spatial logistic regression is performed by the function
  \code{\link{slrm}}. Its result is an object of class \code{"slrm"}.
  There are many methods for this class, including methods for
  \code{print}, \code{fitted}, \code{predict},
  \code{anova}, \code{coef}, \code{logLik}, \code{terms},
  \code{update}, \code{formula} and \code{vcov}. 
  
  For example, if \code{X} is a point pattern (class
  \code{"ppp"}):
  
  \tabular{ll}{
    \emph{command} \tab \emph{model} \cr
    \code{slrm(X ~ 1)} \tab Complete Spatial Randomness \cr
    \code{slrm(X ~ x)} \tab Poisson process with \cr
                \tab intensity loglinear in \eqn{x} coordinate \cr
    \code{slrm(X ~ Z)} \tab Poisson process with \cr
                \tab intensity loglinear in covariate \code{Z}
  }

  \bold{Manipulating a fitted spatial logistic regression}
  
  \tabular{ll}{
    \code{\link{anova.slrm}} \tab Analysis of deviance \cr
    \code{\link{coef.slrm}}  \tab Extract fitted coefficients \cr
    \code{vcov.slrm}  \tab Variance-covariance matrix of fitted coefficients \cr
    \code{\link{fitted.slrm}} \tab Compute fitted probabilities or
    intensity \cr
    \code{\link{logLik.slrm}}   \tab Evalate loglikelihood of fitted
    model \cr
    \code{\link{plot.slrm}}    \tab Plot fitted probabilities or
    intensity \cr
    \code{\link{predict.slrm}} \tab Compute predicted probabilities or
    intensity with new data 
  }
  
  There are many other undocumented methods for this class,
  including methods for \code{print}, \code{update}, \code{formula}
  and \code{terms}. Stepwise model selection is
  possible using \code{step} or \code{stepAIC}.
}


\section{VI. SIMULATION}{

  There are many ways to generate a random point pattern,
  line segment pattern, pixel image or tessellation
  in \pkg{spatstat}. 

  \bold{Random point patterns:}

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
    \code{\link{runifdisc}} \tab
    generate \eqn{n} independent uniform random points in disc\cr
    \code{\link{rstrat}} \tab
    stratified random sample of points \cr
    \code{\link{rsyst}} \tab
    systematic random sample (grid) of points \cr
    \code{\link{rMaternI}}  \tab
    simulate the Mat\'ern Model I inhibition process\cr
    \code{\link{rMaternII}} \tab
    simulate the Mat\'ern Model II inhibition process\cr
    \code{\link{rSSI}} \tab
    simulate Simple Sequential Inhibition process\cr
    \code{\link{rStrauss}} \tab
    simulate Strauss process (perfect simulation)\cr
    \code{\link{rNeymanScott}} \tab
    simulate a general Neyman-Scott process\cr
    \code{\link{rMatClust}} \tab
    simulate the Mat\'ern Cluster process\cr
    \code{\link{rThomas}} \tab
    simulate the Thomas process  \cr
    \code{\link{rGaussPoisson}}  \tab
    simulate the Gauss-Poisson cluster process\cr
    \code{\link{rcell}} \tab
    simulate the Baddeley-Silverman cell process  \cr
    \code{\link{runifpointOnLines}} \tab
    generate \eqn{n} random points along specified line segments \cr
    \code{\link{rpoisppOnLines}} \tab
    generate Poisson random points along specified line segments 
  }

      
  \bold{Resampling a point pattern:}

  \tabular{ll}{
    \code{\link{quadratresample}} \tab block resampling \cr
    \code{\link{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link{rthin}} \tab  random thinning 
  }
  
  \bold{Fitted point process models:}

   If you have fitted a point process model to a point pattern dataset,
   the fitted model can be simulated.

   Cluster process models 
   are fitted by the function \code{\link{kppm}} yielding an
   object of class \code{"kppm"}. To generate one or more simulated
   realisations of this fitted model, use 
   \code{\link{simulate.kppm}}.

   Gibbs point process models 
   are fitted by the function \code{\link{ppm}} yielding an
   object of class \code{"ppm"}. To generate a simulated
   realisation of this fitted model, use \code{\link{rmh}}.
   To generate one or more simulated realisations of the fitted model,
   use \code{\link{simulate.ppm}}.

   \bold{Other random patterns:}

   \tabular{ll}{
     \code{\link{rlinegrid}} \tab
     generate a random array of parallel lines through a window \cr
     \code{\link{rpoisline}} \tab
     simulate the Poisson line process within a window \cr
     \code{\link{rpoislinetess}} \tab
     generate random tessellation using Poisson line process \cr
     \code{\link{rMosaicSet}} \tab
     generate random set by selecting some tiles of a tessellation \cr
     \code{\link{rMosaicField}} \tab
     generate random pixel image by assigning random values
     in each tile of a tessellation
   }

   \bold{Simulation-based inference}

   \tabular{ll}{
    \code{\link{envelope}} \tab critical envelope for Monte Carlo
    test of goodness-of-fit \cr
    \code{\link{qqplot.ppm}} \tab diagnostic plot for interpoint
    interaction
  }
}


\section{VII. TESTS AND DIAGNOSTICS}{

  \bold{Classical hypothesis tests:}
  
  \tabular{ll}{
    \code{\link{quadrat.test}} \tab \eqn{\chi^2}{chi^2} goodness-of-fit
    test on quadrat counts \cr
    \code{\link{kstest}} \tab Kolmogorov-Smirnov goodness-of-fit test\cr
    \code{\link{bermantest}} \tab Berman's goodness-of-fit tests\cr
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


  \bold{Resampling and randomisation procedures}

  You can build your own tests based on randomisation
  and resampling using the following capabilities:
  
  \tabular{ll}{
    \code{\link{quadratresample}} \tab block resampling \cr
    \code{\link{rjitter}} \tab
    apply random displacements to points in a pattern\cr
    \code{\link{rshift}} \tab random shifting of (subsets of) points\cr
    \code{\link{rthin}} \tab  random thinning  
  }
}


\section{VIII. DOCUMENTATION}{
  The online manual entries are quite detailed and should be consulted
  first for information about a particular function.
  
  The paper by Baddeley and Turner (2005a) is a brief overview of
  the package. 
  Baddeley and Turner (2005b) is a more detailed explanation of
  how to fit point process models to data.
  Baddeley (2008) is a complete set of notes from a 2-day workshop
  on the use of \code{spatstat}.

  Type \code{citation("spatstat")} to get these references.
}
\references{
  Baddeley, A. (2008) Analysing spatial point patterns in R.
  Workshop notes. CSIRO online technical publication.
  URL: \code{www.csiro.au/resources/pf16h.html}
  
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
  \email{Adrian.Baddeley@csiro.au}
  \url{http://www.maths.uwa.edu.au/~adrian/}
  and Rolf Turner
  \email{r.turner@auckland.ac.nz}
}
\section{Acknowledgements}{
  Marie-Colette van Lieshout, Rasmus Waagepetersen, 
  Dominic Schuhmacher, Ege Rubak and Kasper Klitgaard Berthelsen
  made substantial contributions of code.
  Additional contributions by
  Ang Qi Wei,
  Sandro Azaele,
  Colin Beale,
  Ricardo Bernhardt,
  Brad Biggerstaff,
  Roger Bivand,
  Florent Bonneu,
  Julian Burgos,
  S. Byers,
  Ya-Mei Chang,
  Jianbao Chen,
  Y.C. Chin,
  Bjarke Christensen,
  Marcelino de la Cruz,
  Peter Dalgaard,
  Peter Diggle,
  Stephen Eglen,
  Agnes Gault,
  Marc Genton,
  Pavel Grabarnik,
  C. Graf,
  Janet Franklin,
  Ute Hahn,
  Andrew Hardegen,
  Mandy Hering, 
  Martin Bogsted Hansen,
  Martin Hazelton,
  Juha Heikkinen,
  Kurt Hornik,
  Ross Ihaka,
  Robert John-Chandran,
  Devin Johnson,
  Jeff Laake,
  Robert Mark,
  Jorge Mateu Mahiques,
  Peter McCullagh,
  Sebastian Wastl Meyer,
  Mi Xiangcheng,
  Jesper Moller,
  Linda Stougaard Nielsen,
  Felipe Nunes,
  Evgeni Parilov,
  Jeff Picka,
  Adrian Raftery,
  Matt Reiter,
  Brian Ripley,
  Barry Rowlingson,
  John Rudge,
  Aila Sarkka,
  Katja Schladitz,
  Bryan Scott,
  Ida-Maria Sintorn,
  Yong Song, 
  Malte Spiess,
  Mark Stevenson,
  Michael Sumner,
  P. Surovy,
  Berwin Turlach,
  Andrew van Burgel,
  Tobias Verbeke,
  Alexendre Villers,
  Hao Wang,
  H. Wendrock,
  Jan Wild
  and
  Selene Wong.
}
\keyword{spatial}
\keyword{package}

