#
#   plot.im.R
#
#  $Revision: 1.27 $   $Date: 2006/08/28 07:11:55 $
#
#  Plotting code for pixel images
#
#  plot.im
#  image.im
#  contour.im
#  persp.im
#
###########################################################################

plot.im <- function(x, ..., ribbon=TRUE, ribsep=0.15, ribwid=0.05, ribn=1024) {
  verifyclass(x, "im")
  main <- deparse(substitute(x))

  imagebreaks <- NULL
  ribbonvalues <- ribbonbreaks <- NULL

  sumry <- summary(x)
  switch(sumry$type,
         real    = {
           vrange <- sumry$range
           trivial <- (diff(vrange) <= .Machine$double.eps)
           ribbonvalues <- seq(vrange[1], vrange[2], length=ribn)
           ribbonrange <- vrange
           ribbonticks <- pretty(ribbonvalues)
           ribbonlabels <- paste(ribbonticks)
         },
         integer = {
           values <- as.vector(x$v)
           values <- values[!is.na(values)]
           uv <- unique(values)
           vrange <- range(uv)
           nvalues <- length(uv)
           trivial <- (nvalues < 2)
           ribbonvalues <- seq(vrange[1], vrange[2], length=ribn)
           ribbonrange <- vrange
           ribbonticks <- pretty(ribbonvalues)
           ribbonlabels <- paste(ribbonticks)
         },
         logical = {
           values <- as.integer(as.vector(x$v))
           values <- values[!is.na(values)]
           trivial <- (length(unique(values)) < 2)
           vrange <- c(0,1)
           imagebreaks <- c(-0.5, 0.5, 1.5)
           ribbonvalues <- c(0,1)
           ribbonrange <- range(imagebreaks)
           ribbonbreaks <- imagebreaks
           ribbonticks <- ribbonvalues
           ribbonlabels <- c("FALSE", "TRUE")
         },
         factor  = {
           lev <- x$lev
           nvalues <- length(lev)
           trivial <- (nvalues < 2)
           # ensure all factor levels plotted separately
           intlev <- as.integer(factor(lev, levels=lev))
           imagebreaks <- c(intlev - 0.5, max(intlev) + 0.5)
           ribbonvalues <- intlev
           ribbonrange <- range(imagebreaks)
           ribbonbreaks <- imagebreaks
           ribbonticks <- ribbonvalues
           ribbonlabels <- paste(lev)
           vrange <- range(intlev)
         },
         stop(paste("Don\'t know how to plot image of type", sQuote(sumry$type)))
         )

  if(!is.null(imagebreaks))
    colourmap <- list(breaks=imagebreaks,
                      col=heat.colors(length(imagebreaks) - 1))
  else
    colourmap <- NULL


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
                 colourmap,
                 list(xlab = "x", ylab = "y"),
                 list(asp = 1, main = main)
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
               colourmap,
               list(xlab = "x", ylab = "y"),
               list(asp = 1, main = main))
  # axes for image
  imax <- resolve.defaults(list(...), list(axes=TRUE))$axes
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
             colourmap)
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
  if(is.null(colmap)) 
    colinfo <- list(col=NULL)
  else {
    # interpret colour map
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
    xbox <- as.rectangle(x)
    zscale <- 0.5 * mean(diff(xbox$xrange), diff(xbox$yrange))/diff(xinfo$range)
  } else
    zscale <- NULL
  
  do.call.matched("persp",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(...),
                                   pop,
                                   colinfo,
                                   list(xlab="x", ylab="y", zlab=xname),
                                   list(scale=FALSE, expand=zscale),
                                   list(main=xname)),
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

