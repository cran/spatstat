#
#	window.S
#
#	A class 'owin' to define the "observation window"
#
#	$Revision: 4.24 $	$Date: 2004/01/08 12:20:05 $
#
#
#	A window may be either
#
#		- rectangular:
#                       a rectangle in R^2
#                       (with sides parallel to the coordinate axes)
#
#		- polygonal:
#			delineated by one or more non-self-intersecting
#                       polygons, possibly including polygonal holes.
#	
#		- digital mask:
#			defined by a binary image
#			whose pixel values are TRUE wherever the pixel
#                       is inside the window
#
#	Any window is an object of class 'owin', 
#       containing at least the following entries:	
#
#		$type:	a string ("rectangle", "polygonal" or  "mask")
#
#		$xrange   
#		$yrange
#			vectors of length 2 giving the real dimensions 
#			of the enclosing box.
#
#	The 'rectangle' type has only these entries.
#
#       The 'polygonal' type has an additional entry
#
#               $bdry
#                       a list of polygons.
#                       Each entry bdry[[i]] determines a closed polygon.
#
#                       bdry[[i]] has components $x and $y which are
#                       the cartesian coordinates of the vertices of
#                       the i-th boundary polygon (without repetition of
#                       the first vertex, i.e. same convention as in the
#                       plotting function polygon().)
#
#
#	The 'mask' type has entries
#
#		$m		logical matrix
#		$dim		its dimension array
#		$xstep,ystep	x and y dimensions of a pixel
#		$xcol	        vector of x values for each column
#               $yrow           vector of y values for each row
#	
#	(the row index corresponds to increasing y coordinate; 
#	 the column index "   "     "   "  "  "  x "   "    ".)
#
#
#-----------------------------------------------------------------------------
#
owin <- function(xrange=c(0,1), yrange=c(0,1), poly=NULL, mask=NULL) {

  ## Exterminate ambiguities
  if(!missing(poly) && !is.null(poly) && !missing(mask) && !is.null(mask))
     stop("Ambiguous -- both polygonal boundary and digital mask supplied")
     
  if(missing(xrange) != missing(yrange))
    stop("If one of xrange, yrange is specified then both must be.")

  if(missing(poly) && missing(mask)) {
    ######### rectangle #################
    if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] <= xrange[1])
      stop("xrange should be a vector of length 2 giving (xmin, xmax)")
    if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] <= yrange[1])
      stop("yrange should be a vector of length 2 giving (ymin, ymax)")
    w <- list(type="rectangle", xrange=xrange, yrange=yrange)
    class(w) <- "owin"
    return(w)
  } else if(!missing(poly)) {
    ######### polygonal boundary ########
    #
    # test whether it's a single polygon or multiple polygons
    if(verify.xypolygon(poly, fatal=FALSE))
      psingle <- TRUE
    else if(all(unlist(lapply(poly, verify.xypolygon, fatal=FALSE))))
      psingle <- FALSE
    else
      stop("poly must be either a list(x,y) or a list of list(x,y)")
                  
    if(psingle) {
      # single boundary polygon
      if(area.xypolygon(poly) < 0)
        stop("Area of polygon is negative - maybe traversed in wrong direction?")
      bdry <- list(poly)
    } else {
      # multiple boundary polygons
      bdry <- poly
      if(sum(unlist(lapply(poly, area.xypolygon))) < 0)
        stop(paste("Area of window is negative;\n",
             "check that all polygons were traversed in the right direction"))
    }

    actual.xrange <- range(unlist(lapply(bdry, function(a) a$x)))
    if(missing(xrange))
      xrange <- actual.xrange
    else {
      if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] <= xrange[1])
        stop("xrange should be a vector of length 2 giving (xmin, xmax)")
      if(!all(xrange == range(c(xrange, actual.xrange))))
        stop("polygon's x coordinates outside xrange")
    }
    
    actual.yrange <- range(unlist(lapply(bdry, function(a) a$y)))
    if(missing(yrange))
      yrange <- actual.yrange
    else {
      if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] <= yrange[1])
        stop("yrange should be a vector of length 2 giving (ymin, ymax)")
      if(!all(yrange == range(c(yrange, actual.yrange))))
      stop("polygon's y coordinates outside yrange")
    }

    w <- list(type="polygonal", xrange=xrange, yrange=yrange, bdry=bdry)
    class(w) <- "owin"
    return(w)
    
  } else if(!missing(mask)) {
    ######### digital mask #####################
    
    if(!is.matrix(mask))
      stop("\`mask\' must be a matrix")
    if(!is.logical(mask))
      stop("The entries of \`mask\' must be logical")
    
    nc <- ncol(mask)
    nr <- nrow(mask)

    if(missing(xrange) && missing(yrange)) {
      # take pixels to be 1 x 1 unit
      xrange <- c(0,nc)
      yrange <- c(0,nr)
    } else {
      if(!is.vector(xrange) || length(xrange) != 2 || xrange[2] <= xrange[1])
        stop("xrange should be a vector of length 2 giving (xmin, xmax)")
      if(!is.vector(yrange) || length(yrange) != 2 || yrange[2] <= yrange[1])
        stop("yrange should be a vector of length 2 giving (ymin, ymax)")
    }

    xstep <- diff(xrange)/nc
    ystep <- diff(yrange)/nr
  
    out <- list(type     = "mask",
                xrange   = xrange,
                yrange   = yrange,
                dim      = c(nr, nc),
                xstep    = xstep,
                ystep    = ystep,
                warnings = c(
"Row index corresponds to increasing y coordinate; column to increasing x",
"Transpose matrices to get the standard presentation in S",
"Example: image(result$xcol,result$yrow,t(result$d))"
                ),
                xcol    = seq(xrange[1]+xstep/2, xrange[2]-xstep/2, length=nc),
                yrow    = seq(yrange[1]+ystep/2, yrange[2]-ystep/2, length=nr),
                m       = mask)
    class(out) <- "owin"
    return(out)
  }
  # never reached
  NULL
}

#
#-----------------------------------------------------------------------------
#

is.owin <- function(x) { inherits(x, "owin") }

as.owin <- function(W) {
	# Tries to interpret data as an object of class 'window'
	# W may be
	#	an object of class 'window'
	#	a structure with entries xrange, yrange
	#	a four-element vector (interpreted xmin, xmax, ymin, ymax)
	#	a structure with entries xl, xu, yl, yu
	#	an object of class 'ppp'
	#	an object of class 'im'

	if(verifyclass(W, "owin", fatal=FALSE))
		return(W)
        else if(verifyclass(W, "ppp", fatal=FALSE))
		return(W$window)
        else if(verifyclass(W, "im", fatal=FALSE)) 
                return(owin(W$xrange, W$yrange, mask=!is.na(W$v)))
        else if(checkfields(W, c("xrange", "yrange"))) {
		Z <- owin(W$xrange, W$yrange)
		return(Z)
	} else if(is.vector(W) && is.numeric(W) && length(W) == 4) {
		Z <- owin(W[1:2], W[3:4])
		return(Z)
	} else if(checkfields(W, c("xl", "xu", "yl", "yu"))) {
                W <- as.list(W)
		Z <- owin(c(W$xl, W$xu),c(W$yl, W$yu))
		return(Z)
        } else if(checkfields(W, c("x", "y", "area"))
                  && checkfields(W$area, c("xl", "xu", "yl", "yu"))) {
                V <- as.list(W$area)
                Z <- owin(c(V$xl, V$xu),c(V$yl, V$yu))
                return(Z)
	} else
		stop("Can't interpret W as a window")
}		

#
#-----------------------------------------------------------------------------
#
#
as.rectangle <- function(...) {
        w <- as.owin(...)
        return(owin(w$xrange, w$yrange))
}

#
#-----------------------------------------------------------------------------
#
as.mask <- function(w, eps=NULL, dimyx=NULL) {
#	eps:		   grid mesh (pixel) size
#	dimyx:		   dimensions of pixel raster
#
#  NOTE: rows <=> y coordinate

        verifyclass(w, "owin")
        if(w$type == "mask")
          return(w)

        # determine row & column dimensions of output raster
        if(!is.null(dimyx)) {
          nr <- dimyx[1]
          nc <- dimyx[2]
        } else {
          # use pixel size 'eps'
          if(!is.null(eps)) {
            nr <- ceiling(diff(w$yrange)/eps)
            nc <- ceiling(diff(w$xrange)/eps)
          } else {
            # use spatstat defaults
            np <- spatstat.options("npixel")[[1]]
            if(!is.numeric(np) || length(np) > 2)
                   stop("Illegal value for spatstat.options(\"npixel\")")
            if(length(np) == 1)
                     nr <- nc <- np[1]
            else {
                     nr <- np[2]  
                     nc <- np[1]
            }
          }
        }

        # mask with all entries True
        out <- owin(w$xrange, w$yrange,
                    mask=matrix(TRUE, nrow=nr, ncol=nc))
        
        switch(w$type,
               rectangle={
                 return(out)
               },
               polygonal={
                 # test every pixel
                 x <- as.vector(raster.x(out))
                 y <- as.vector(raster.y(out))
                 value <- inside.owin(x, y, w)
                 return(owin(w$xrange, w$yrange,
                             mask=matrix(value, nrow=nr, ncol=nc)))
               },
               stop("Unrecognised window type \`", w$type, "\'")
               )
}

#
#-----------------------------------------------------------------------------
#
as.polygonal <- function(W) {
  # the only sensible use of this function is to
  # convert a rectangle to a polygon
  verifyclass(W, "owin")
  switch(W$type,
         rectangle = {
           xr <- W$xrange
           yr <- W$yrange
           return(owin(xr, yr, poly=list(x=xr[c(1,2,2,1)],y=yr[c(1,1,2,2)])))
         },
         polygonal = {
           return(W)
         },
         mask = {
           stop("A mask cannot be converted to a polygon")
         }
         )
}

#
# ----------------------------------------------------------------------

validate.mask <- function(w, fatal=TRUE) {
  verifyclass(w, "owin", fatal=fatal)
  if(w$type == "mask")
    return(TRUE)
  if(fatal)
      stop(deparse(substitute(w)), "is not a binary mask")
  else {
      warning(deparse(substitute(w)), "is not a binary mask")
      return(FALSE)
  }
}
             
raster.x <- function(w) {
	validate.mask(w)
        di <- w$dim
        nr <- di[1]
        nc <- di[2]
	m <- matrix(0, nrow=nr, ncol=nc)
	matrix(w$xcol[col(m)], nrow=nr, ncol=nc)
}

raster.y <- function(w) {
	validate.mask(w)
        di <- w$dim
        nr <- di[1]
        nc <- di[2]
	m <- matrix(0, nrow=nr, ncol=nc)
	matrix(w$yrow[row(m)], nrow=nr, ncol=nc)
}

nearest.raster.point <- function(x,y,w) {
	validate.mask(w)
	nr <- w$dim[1]
	nc <- w$dim[2]
	cc <- round(0.5 + (x - w$xrange[1])/w$xstep)
	rr <- round(0.5 + (y - w$yrange[1])/w$ystep)
	cc <- pmax(1,pmin(cc, nc))
	rr <- pmax(1,pmin(rr, nr))
	return(list(row=rr, col=cc))
}

#------------------------------------------------------------------
		
bounding.box <- function(...) {
        wins <- list(...)

        if(length(wins) == 0)
          stop("No arguments supplied")

        if(length(wins) > 1) {
          # multiple arguments -- compute bounding box for each
          boxes <- lapply(wins, bounding.box)
          # take bounding box of these boxes
          xrange <- range(unlist(lapply(boxes, function(b){b$xrange})))
          yrange <- range(unlist(lapply(boxes, function(b){b$yrange})))
          return(owin(xrange, yrange))
        }

        w <- wins[[1]]
          
        # determine a tight bounding box for the window w
        verifyclass(w, "owin")

        switch(w$type,
               rectangle = {
                 return(w)
               },
               polygonal = {
                 bdry <- w$bdry
                 xr <- range(unlist(lapply(bdry, function(a) a$x)))
                 yr <- range(unlist(lapply(bdry, function(a) a$y)))
                 return(owin(xr, yr))
               },
               mask = {
                 m <- w$m
                 x <- raster.x(w)
                 y <- raster.y(w)
                 xr <- range(x[m])
                 yr <- range(y[m])
                 return(owin(xr, yr))
               },
               stop("unrecognised window type", w$type)
               )
}
  
complement.owin <- function(w) {
	verifyclass(w, "owin")
        switch(w$type,
               mask = {
                 w$m <- !(w$m)
               },
               polygonal = {

                 bdry <- w$bdry
                 
                 # bounding box, in anticlockwise order
                 box <- list(x=w$xrange[c(1,2,2,1)],
                             y=w$yrange[c(1,1,2,2)])
                 boxarea <- area.xypolygon(box)
                 
                 # first check whether one of the current boundary polygons
                 # is the bounding box itself (with + sign)
                 nvert <- unlist(lapply(bdry, function(a) { length(a$x) }))
                 area <- unlist(lapply(bdry, area.xypolygon))
                 boxarea.mineps <- boxarea * (1 - .Machine$single.eps)
                 is.box <- (nvert == 4 & area >= boxarea.mineps)
                 if(sum(is.box) > 1)
                   stop("Internal error: multiple copies of bounding box")
                 
                 # if box is present (with + sign), remove it
                 if(any(is.box))
                   bdry <- bdry[!is.box]
                 
                 # reverse the direction of each polygon
                 bdry <- lapply(bdry, reverse.xypolygon)
                 
                 # if box was absent, add it
                 if(!any(is.box))
                   bdry <- c(bdry, list(box))   # sic
                 
                 # put back into w
                 w$bdry <- bdry
               },
               rectangle = {
                 stop("window is a rectangle - its complement is empty")
               },
               stop("unrecognised window type", w$type)
               )
	return(w)
}

#-----------------------------------------------------------

inside.owin <- function(x, y, w) {
  # test whether (x,y) is inside window w
  # x, y may be vectors 
  
  verifyclass(w, "owin")

  # test whether inside bounding rectangle
  xr <- w$xrange
  yr <- w$yrange
  frameok <- (xr[1] <= x) & (x <= xr[2]) & (yr[1] <= y) & (y <= yr[2])

  if(all(!frameok))  # all points OUTSIDE window - no further work needed
    return(frameok)

  ok <- frameok
  switch(w$type,
         rectangle = {
           return(ok)
         },
         polygonal = {
           xy <- list(x=x,y=y)
           bdry <- w$bdry
           total <- rep(0, length(x))
           on.bdry <- rep(FALSE, length(x))
           for(i in seq(bdry)) {
             score <- inside.xypolygon(xy, bdry[[i]], test01=FALSE)
             total <- total + score
             on.bdry <- on.bdry | attr(score, "on.boundary")
           }
           # any points identified as belonging to the boundary get score 1
           total[on.bdry] <- 1
           # check for sanity now..
           if(any(total * (1-total) != 0)) 
             stop("internal error: some total scores are neither 0 nor 1")
           return(ok & (total != 0))
         },
         mask = {
           # consider only those points which are inside the frame
           xf <- x[frameok]
           yf <- y[frameok]
           # map locations to raster (row,col) coordinates
           loc <- nearest.raster.point(xf,yf,w)
           # look up mask values
           mas <- w$m
           nf <- sum(frameok)
           okf <- logical(nf)
           for(i in 1:nf) 
             okf[i] <- mas[loc$row[i],loc$col[i]]
           # insert into 'ok' vector
           ok[frameok] <- okf
           return(ok)
         },
         stop("unrecognised window type", w$type)
         )
}

#-------------------------------------------------------------------------
  
print.owin <- function(x, ...) {
  verifyclass(x, "owin")
  cat("window: ")
  switch(x$type,
         rectangle={
           cat("rectangle = ")
         },
         polygonal={
           cat("polygonal boundary\n")
           cat("enclosing rectangle: ")
         },
         mask={
           cat("binary image mask\n")
           cat("enclosing rectangle: ")
         }
         )
    cat(paste("[",
              x$xrange[1],
              ",",
              x$xrange[2],
              "] x [",
              x$yrange[1],
              ",",
              x$yrange[2],
              "]\n"))
}

summary.owin <- function(object, ...) {
  verifyclass(object, "owin")
  result <- list(xrange=object$xrange,
                 yrange=object$yrange,
                 type=object$type,
                 area=area.owin(object))
  switch(object$type,
         rectangle={
         },
         polygonal={
           poly <- object$bdry
           result$npoly <- npoly <- length(poly)
           if(npoly == 1) {
             result$areas <- area.xypolygon(poly[[1]])
             result$nvertices <- length(poly[[1]]$x)
           } else {
             result$areas <- unlist(lapply(poly, area.xypolygon))
             result$nvertices <- unlist(lapply(poly,
                                               function(a) {length(a$x)}))
           }
           result$nhole <- sum(result$areas < 0)
         },
         mask={
           result$npixels <- object$dim
         }
         )
  class(result) <- "summary.owin"
  result
}

print.summary.owin <- function(x, ...) {
  verifyclass(x, "summary.owin")
  cat("Window: ")
  switch(x$type,
         rectangle={
           cat("rectangle = ")
         },
         polygonal={
           cat("polygonal boundary\n")
           if(x$npoly == 1) {
             cat(paste("single connected closed polygon with",
                       x$nvertices, 
                       "vertices\n"))
           } else {
             cat(paste(x$npoly, "separate polygons ("))
             if(x$nhole == 0) cat("no holes)\n")
             else if(x$nhole == 1) cat("1 hole)\n")
             else cat(paste(x$nhole, "holes)\n"))
             print(data.frame(vertices=x$nvertices,
                              area=x$areas,
                              relative.area=signif(x$areas/x$area,3),
                              row.names=paste("polygon",
                                1:(x$npoly),
                                ifelse(x$areas < 0, "(hole)", "")
                                )))
           }
           cat("enclosing rectangle: ")
         },
         mask={
           cat("binary image mask\n")
           di <- x$npixels
           cat(paste(di[1], "x", di[2], "pixel array\n"))
           cat("enclosing rectangle: ")
         }
         )
  cat(paste("[",
            x$xrange[1],
            ",",
            x$xrange[2],
            "] x [",
            x$yrange[1],
            ",",
            x$yrange[2],
            "]\n"))
  cat(paste("Window area = ", x$area, "\n"))
  return(invisible(x))
}

  
