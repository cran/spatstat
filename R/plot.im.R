#
#   plot.im.R
#
#  $Revision: 1.3 $   $Date: 2005/02/08 18:10:22 $
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
  main <- deparse(substitute(x))
  verifyclass(x, "im")
  vrange <- summary(x)$range
  add <- resolve.defaults(list(...), list(add=FALSE))$add
  if(!identical(ribbon, TRUE)
     || identical(add, TRUE)
     || diff(vrange) <= .Machine$double.eps)
    {
      # plot image without ribbon
      do.call("image",
               resolve.defaults(list(x$xcol, x$yrow, t(x$v)),
               list(...),
               list(xlab = "x", ylab = "y", asp = 1, main = main)))
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
  do.call("plot",
          resolve.defaults(list(0, 0,  type="n", axes=FALSE, asp=1,
                            xlim=bb.all$xrange, ylim=bb.all$yrange),
                           list(...), list(main=main, xlab="", ylab="")))
  # plot image
    do.call("image", resolve.defaults(list(x$xcol, x$yrow, t(x$v), add=TRUE),
        list(...), list(xlab = "x", ylab = "y", asp = 1, main = main)))
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
  values <- seq(vrange[1], vrange[2], length=ribn)
  pv <- pretty(values)
  values <- matrix(values, nrow=1)
  ycoords <- seq(bb.rib$yrange[1], bb.rib$yrange[2], length=ribn+1)
  do.call("image",
          resolve.defaults(
              list(bb.rib$xrange, ycoords, values, add=TRUE),
              list(...)))
  plot(as.owin(bb.rib), add=TRUE)
  # ticks for ribbon image
  scal <- diff(bb.rib$yrange)/diff(vrange)
  at.y <- bb.rib$yrange[1] + scal * (pv - vrange[1])
  par(yaxp=c(bb.rib$yrange, length(pv)))
  axis(4, at=at.y, labels=paste(pv), pos=bb.rib$xrange[2])
  # return image call
  return(invisible(NULL))
} 

########################################################################

image.im <- plot.im

########################################################################

persp.im <- function(x, ...) {
  xname <- deparse(substitute(x))
  do.call("persp",
          resolve.defaults(list(x$xcol, x$yrow, t(x$v)),
                           list(...),
                           list(xlab="x", ylab="y", zlab=xname),
                           list(main=xname)))
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
