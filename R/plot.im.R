#
#   plot.im.R
#
#  $Revision: 1.6 $   $Date: 2005/12/05 08:40:40 $
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

  if(!identical(ribbon, TRUE)
     || identical(add, TRUE)
     || trivial)
    {
      # plot image without ribbon
      do.call.matched("image.default",
               resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                list(...),
                                colourmap,
                                list(xlab = "x", ylab = "y"),
                                list(asp = 1, main = main)))
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
    do.call.matched("image.default",
                    resolve.defaults(
                                     list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                     list(add=TRUE),
                                     list(...),
                                     colourmap,
                                     list(xlab = "x", ylab = "y"),
                                     list(asp = 1, main = main)))
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
  do.call.matched("image.default",
                  resolve.defaults(
                                   list(x=bb.rib$xrange, y=ycoords,
                                        z=matrix(ribbonvalues, nrow=1),
                                        add=TRUE),
                                   list(...),
                                   colourmap))
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

persp.im <- function(x, ...) {
  xname <- deparse(substitute(x))
  if(summary(x)$type == "factor")
    stop("Perspective plot is inappropriate for factor-valued image")
  do.call.matched("persp",
                  resolve.defaults(list(x=x$xcol, y=x$yrow, z=t(x$v)),
                                   list(...),
                                   list(xlab="x", ylab="y", zlab=xname),
                                   list(main=xname)),
                  funargs=list("x", "y", "z",
                               "xlim", "ylim", "zlim",
                               "xlab", "ylab", "zlab",
                               "main", "sub",
                               "theta", "phi", "r", "d", "scale",
                               "expand", "col", "border",
                               "ltheta", "lphi", "shade", "box",
                               "axes", "nticks", "ticktype"))
}


######################################################################

contour.im <- function (x, ...)
{
  main <- deparse(substitute(x))
  add <- resolve.defaults(list(...), list(add=FALSE))$add
  if(!add) 
    do.call("plot",
            resolve.defaults(list(range(x$xcol), range(x$yrow), type="n"),
                             list(...),
                             list(asp = 1, xlab="x", ylab="y", main=main)))
  do.call("contour",
          resolve.defaults(list(x$xcol, x$yrow, t(x$v), add=TRUE),
                           list(...)))
}
