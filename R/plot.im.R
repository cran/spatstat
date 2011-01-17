#
#   plot.im.R
#
#  $Revision: 1.47 $   $Date: 2010/12/13 04:47:50 $
#
#  Plotting code for pixel images
#
#  plot.im
#  image.im
#  contour.im
#  persp.im
#
###########################################################################

plot.im <- function(x, ...,
                    col=NULL,
                    ribbon=TRUE, ribsep=0.15, ribwid=0.05, ribn=1024,
                    ribscale=1) {
  verifyclass(x, "im")
  main <- deparse(substitute(x))

  zlim <- list(...)$zlim

  imagebreaks <- NULL
  ribbonvalues <- ribbonbreaks <- NULL

  clamp <- function(x, v, tol=0.02 * diff(v)) {
    ok <- (x >= v[1] - tol) & (x <= v[2] + tol)
    x[ok]
  }

  colmap <- if(inherits(col, "colourmap")) col else NULL
  
  sumry <- summary(x)

  switch(sumry$type,
         real    = {
           vrange <- sumry$range
           vrange <- range(zlim, vrange)
           if(!is.null(colmap)) {
             # explicit colour map
             s <- summary(colmap)
             if(s$discrete)
               stop("Discrete colour map is not applicable to real values")
             imagebreaks <- s$breaks
             vrange <- range(imagebreaks)
             col <- s$outputs
           }
           trivial <- (diff(vrange) <= .Machine$double.eps)
           if(!trivial) {
             ribbonvalues <- seq(vrange[1], vrange[2], length=ribn)
             ribbonrange <- vrange
             ribbonticks <- clamp(pretty(ribscale * ribbonvalues)/ribscale,
                                  vrange)
             ribbonlabels <- paste(ribbonticks * ribscale)
           }
         },
         integer = {
           values <- as.vector(x$v)
           values <- values[!is.na(values)]
           uv <- unique(values)
           vrange <- range(uv)
           vrange <- range(zlim, vrange)
           nvalues <- length(uv)
           trivial <- (nvalues < 2)
           if(!trivial){
             ribbonticks <- clamp(pretty(vrange * ribscale)/ribscale, vrange)
             ribbonticks <- ribbonticks[ribbonticks %% 1 == 0]
             if(identical(all.equal(ribbonticks,
                                    vrange[1]:vrange[2]), TRUE)) {
               # each possible value will be printed on ribbon label
               ribbonvalues <- vrange[1]:vrange[2]
               imagebreaks <- c(ribbonvalues - 0.5, vrange[2] + 0.5)
               ribbonrange <- range(imagebreaks)
               ribbonticks <- ribbonvalues
               ribbonlabels <- paste(ribbonticks * ribscale)
             } else {
               # not all values will be printed on ribbon label
               ribn <- min(ribn, diff(vrange)+1)
               ribbonvalues <- seq(vrange[1], vrange[2], length=ribn)
               ribbonrange <- vrange
               ribbonlabels <- paste(ribbonticks * ribscale)
             }
           }
           if(!is.null(colmap)) {
             # explicit colour map
             s <- summary(colmap)
             if(!s$discrete)
               imagebreaks <- s$breaks
             else
               imagebreaks <- c(s$inputs[1] - 0.5, s$inputs + 0.5)
             col <- s$outputs
           }
         },
         logical = {
           values <- as.integer(as.vector(x$v))
           values <- values[!is.na(values)]
           uv <- unique(values)
           trivial <- (length(uv) < 2)
           vrange <- c(0,1)
           imagebreaks <- c(-0.5, 0.5, 1.5)
           ribbonvalues <- c(0,1)
           ribbonrange <- range(imagebreaks)
           ribbonbreaks <- imagebreaks
           ribbonticks <- ribbonvalues
           ribbonlabels <- c("FALSE", "TRUE")
           if(!is.null(colmap)) 
             col <- colmap(c(FALSE,TRUE))
         },
         factor  = {
           lev <- levels(x)
           nvalues <- length(lev)
           trivial <- (nvalues < 2)
           # ensure all factor levels plotted separately
           fac <- factor(lev, levels=lev)
           intlev <- as.integer(fac)
           imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
           ribbonvalues <- intlev
           ribbonrange <- range(imagebreaks)
           ribbonbreaks <- imagebreaks
           ribbonticks <- ribbonvalues
           ribbonlabels <- paste(lev)
           vrange <- range(intlev)
           if(!is.null(colmap)) 
             col <- colmap(fac)
         },
         stop(paste("Do not know how to plot image of type", sQuote(sumry$type)))
         )

  
  # determine colour map
  if(!is.null(colmap)) {
    # explicit colour map object
    colourinfo <- list(breaks=imagebreaks, col=col)
  } else {
    # compile colour information using default colour values
    colfun <- spatstat.options("image.colfun")
    if(!is.null(imagebreaks)) 
      colourinfo <- list(breaks=imagebreaks, col=colfun(length(imagebreaks)-1))
    else 
      colourinfo <- list(col=colfun(255))
    # user-specified colour values
    if(!is.null(col)) {
      colourinfo$col <- col
      if(!is.null(colourinfo$breaks)) {
      # check consistency
        nvalues <- length(colourinfo$breaks) - 1
        if(length(col) != nvalues)
          stop(paste("Length of argument", dQuote("col"),
                     paren(paste(length(col))),
                     "does not match the number of distinct values",
                     paren(paste(nvalues))))
      }
    }
  }
  
  add <- resolve.defaults(list(...), list(add=FALSE))$add

  image.doit <- function(...) {
    do.call.matched("image.default",
                    resolve.defaults(...),
                    extrargs=c("main", "asp", "sub", "axes", "ann"))
  }
  
  if(!identical(ribbon, TRUE)
     || identical(add, TRUE)
     || trivial)
    {
      # plot image without ribbon
      image.doit(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                 list(...),
                 colourinfo,
                 list(xlab = "", ylab = ""),
                 list(asp = 1, main = main, axes=FALSE)
                 )
      return(invisible(NULL))
    }
  # determine plot region
  # image at left, ribbon at right
  bb <- owin(x$xrange, x$yrange)
  xwidth <- diff(bb$xrange)
  xheight <- diff(bb$yrange)
  xsize <- max(xwidth, xheight)
  bb.rib <- owin(bb$xrange[2] + c(ribsep, ribsep+ribwid) * xsize,
                 bb$yrange)
  bb.all <- bounding.box(bb.rib, bb)
  # establish coordinate system
  do.call.matched("plot.default",
          resolve.defaults(list(x=0, y=0,  type="n", axes=FALSE, asp=1,
                            xlim=bb.all$xrange, ylim=bb.all$yrange),
                           list(...), list(main=main, xlab="", ylab="")))
  # plot image
    image.doit(list(x=x$xcol, y=x$yrow, z=t(x$v)),
               list(add=TRUE),
               list(...),
               colourinfo,
               list(xlab = "", ylab = ""),
               list(asp = 1, main = main))
  # axes for image
  imax <- resolve.defaults(list(...), list(axes=FALSE))$axes
  if(imax) {
    px <- pretty(bb$xrange)
    py <- pretty(bb$yrange)
    axis(1, at=px, pos=bb$yrange[1])
    axis(2, at=py, pos=bb$xrange[1])
    rect(x$xrange[1], x$yrange[1], x$xrange[2], x$yrange[2])
  }
  # plot ribbon image containing the range of image values
  ycoords <- seq(bb.rib$yrange[1], bb.rib$yrange[2],
                 length=length(ribbonvalues)+1)
  image.doit(list(x=bb.rib$xrange, y=ycoords,
                  z=matrix(ribbonvalues, nrow=1),
                  add=TRUE, main="", sub=""),
             list(...),
             colourinfo)
  plot(as.owin(bb.rib), add=TRUE)
  # ticks for ribbon image
  scal <- diff(bb.rib$yrange)/diff(ribbonrange)
  at.y <- bb.rib$yrange[1] + scal * (ribbonticks - ribbonrange[1])
  par(yaxp=c(bb.rib$yrange, length(ribbonticks)))
  axis(4, at=at.y, labels=ribbonlabels, pos=bb.rib$xrange[2])
  #
  return(invisible(NULL))
} 


########################################################################

image.im <- plot.im

########################################################################

persp.im <- function(x, ..., colmap=NULL) {
  xname <- deparse(substitute(x))
  xinfo <- summary(x)
  if(xinfo$type == "factor")
    stop("Perspective plot is inappropriate for factor-valued image")
  pop <- spatstat.options("par.persp")
  # check for common error
  if(!is.null(col <- list(...)$col) && !is.matrix(col))
    warning("Argument col is not a matrix. Did you mean colmap?")
  # colour map?
  if(is.null(colmap)) {
    colinfo <- list(col=NULL)
  } else if(inherits(colmap, "colourmap")) {
    # colour map object
    # apply colour function to image data
    colval <- eval.im(colmap(x))
    colval <- t(as.matrix(colval))
    # strip one row and column for input to persp.default
    colval <- colval[-1, -1]
    # replace NA by arbitrary value
    isna <- is.na(colval)
    if(any(isna)) {
      stuff <- attr(colmap, "stuff")
      colvalues <- stuff$outputs
      colval[isna] <- colvalues[1]
    }
    # pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  } else {
    # interpret 'colmap' as colour map
    if(is.list(colmap) && all(c("breaks", "col") %in% names(colmap))) {
      breaks <- colmap$breaks
      colvalues <- colmap$col
    } else if(is.vector(colmap)) {
      colvalues <- colmap
      breaks <- quantile(x, seq(0,1,length=length(colvalues)+1))
      if(!all(ok <- !duplicated(breaks))) {
        breaks <- breaks[ok]
        colvalues <- colvalues[ok[-1]]
      }
    } else warning("Unrecognised format for colour map")
    # apply colour map to image values
    colid <- cut.im(x, breaks=breaks, include.lowest=TRUE)
    colval <- eval.im(colvalues[unclass(colid)])
    colval <- t(as.matrix(colval))
    nr <- nrow(colval)
    nc <- ncol(colval)
    colval <- colval[-1, -1]
    colval[is.na(colval)] <- colvalues[1]
    # pass colour matrix (and suppress lines)
    colinfo <- list(col=colval, border=NA)
  }

  # get reasonable z scale while fixing x:y aspect ratio
  if(xinfo$type %in% c("integer", "real")) {
    zrange <- xinfo$range
    if(diff(zrange) > 0) {
      xbox <- as.rectangle(x)
      zscale <- 0.5 * mean(diff(xbox$xrange), diff(xbox$yrange))/diff(zrange)
      zlim <- zrange
    } else {
      zscale <- NULL
      zlim <- c(0,2) * xinfo$mean
    }
  } else 
    zscale <- zlim <- NULL

  do.call.matched("persp",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(...),
                                   pop,
                                   colinfo,
                                   list(xlab="x", ylab="y", zlab=xname),
                                   list(scale=FALSE, expand=zscale, zlim=zlim),
                                   list(main=xname),
                                   .StripNull=TRUE),
                  funargs=.Spatstat.persp.args)
}

.Spatstat.persp.args <- list("x", "y", "z",
                             "xlim", "ylim", "zlim",
                             "xlab", "ylab", "zlab",
                             "main", "sub",
                             "theta", "phi", "r", "d", "scale",
                             "expand", "col", "border",
                             "ltheta", "lphi", "shade", "box",
                             "axes", "nticks", "ticktype")

######################################################################

contour.im <- function (x, ..., main, axes=TRUE, add=FALSE)
{
  defaultmain <- deparse(substitute(x))
  sop <- spatstat.options("par.contour")
  if(missing(main)) 
    main <- resolve.defaults(sop, list(main=defaultmain))$main
  if(missing(add))
    add <- resolve.defaults(sop, list(add=FALSE))$add
  if(missing(axes))
     axes <- resolve.defaults(sop, list(axes=TRUE))$axes
  if(!add) {
    # new plot
    if(axes) # with axes
      do.call.matched("plot.default",
                      resolve.defaults(
                                       list(x = range(x$xcol),
                                            y = range(x$yrow),
                                            type = "n"),
                                       list(...),
                                       list(asp = 1, xlab = "x",
                                            ylab = "y", main = main)))
    else { # box without axes
      rec <- owin(x$xrange, x$yrange)
      do.call.matched("plot.owin",
                      resolve.defaults(list(x=rec),
                                       list(...),
                                       list(main=main)))
    }
  }
  do.call.matched("contour.default",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(add=TRUE),
                                   list(...)))
}

