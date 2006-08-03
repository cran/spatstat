#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.16 $	$Date: 2006/08/02 06:35:19 $
#
#
#

plot.owin <- function(x, main, add=FALSE, ..., box=TRUE, edge=0.04)
{
#
# Function plot.owin.  A method for plot.
#
  if(missing(main))
    main <- deparse(substitute(x))
  W <- x
  verifyclass(W, "owin")

####
  main <- paste(main)
  if(length(main) > 1)
    main <- paste(main, collapse="\n")
  if(nchar(main) == 0)
    nlines <- 0
  else
    nlines <- length(strsplit(main, "\n")[[1]])

#########        
  xlim <- xr <- W$xrange
  ylim <- yr <- W$yrange

####################################################  
  if(!add) {
    # new plot
    # allow space for main title
    guesslinespace <- 0.08 * diff(yr)
    ylim[2] <- ylim[2] + nlines * guesslinespace
    # set up plot with equal scales
    do.call.matched("plot.default",
                    resolve.defaults(list(x=numeric(0), y=numeric(0),
                                          type="n"),
                                     list(...),
                                     list(xlim=xlim, ylim=ylim,
                                          ann=FALSE, axes=FALSE,
                                          asp=1.0)))
    # add title in a reasonable place!
    if(nlines > 0) {
      parval <- resolve.defaults(list(...), par())
      mainheight <- strheight(main, units="user", cex=parval$cex.main)
      gapheight <- (strheight("b\nb", units="user", cex=parval$cex.main)
                    - 2 * strheight("b", units="user", cex=parval$cex.main))
      text(x=mean(xr), y=yr[2] + mainheight + 0.5 * gapheight, labels=main,
           cex=parval$cex.main,
           col=parval$col.main,
           font=parval$font.main)
    }
  }
# Draw window

  switch(W$type,
         rectangle = {
         },
         polygonal = {
           p <- W$bdry
           for(i in seq(p))
             do.call.matched("polygon",
                             resolve.defaults(
                                              list(x=p[[i]]),
                                              list(...),
                                              list(density=0)),
                             extrargs="lwd")
         },
         mask = {
           do.call.matched("image.default",
                           resolve.defaults(
                           list(x=W$xcol, y=W$yrow, z=!t(W$m), add=TRUE),
                           list(...),
                           spatstat.options("par.binary")))
         },
         stop(paste("Don't know how to plot window of type", sQuote(W$type)))
         )

# Draw surrounding box
  wantbox <- (!missing(box) && box) || (missing(box) && W$type == "rectangle")
  if(wantbox)
    do.call.matched("segments",
                    resolve.defaults(
                                     list(x0=xr[c(1,2,2,1)],
                                          y0=yr[c(1,1,2,2)],
                                          x1=xr[c(2,2,1,1)],
                                          y1=yr[c(1,2,2,1)]),
                                     list(...)))
  invisible()
}





