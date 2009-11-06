#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.29 $	$Date: 2009/11/06 00:40:32 $
#
#
#

plot.owin <- function(x, main, add=FALSE, ..., box, edge=0.04,
                      hatch=FALSE, angle=45, spacing=diameter(x)/50)
{
#
# Function plot.owin.  A method for plot.
#
  if(missing(main))
    main <- deparse(substitute(x))
  W <- x
  verifyclass(W, "owin")
  if(missing(box))
    box <- (W$type == "mask")
  else
    stopifnot(is.logical(box) && length(box) == 1)

####
  if(is.expression(main)) 
    nlines <- 1
  else {
    # convert to string and count number of lines
    main <- paste(main)
    if(length(main) > 1)
      main <- paste(main, collapse="\n")
    if(nchar(main) == 0)
      nlines <- 0
    else
      nlines <- length(strsplit(main, "\n")[[1]])
  }
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
  
# Draw surrounding box
  if(box)
    do.call.matched("segments",
                    resolve.defaults(
                                     list(x0=xr[c(1,2,2,1)],
                                          y0=yr[c(1,1,2,2)],
                                          x1=xr[c(2,2,1,1)],
                                          y1=yr[c(1,2,2,1)]),
                                     list(...)))
  
# If type = "n", do not plot the window.
    type <- resolve.defaults(list(...), list(type=NULL))$type
    if(!is.null(type) && type == "n")
      return(invisible(NULL))

  
# Draw window

  switch(W$type,
         rectangle = {
           Wpoly <- as.polygonal(W)
           po <- Wpoly$bdry[[1]]
           do.call.matched("polygon",
                           resolve.defaults(list(x=po),
                                            list(...)),
                           extrargs="lwd")
           if(hatch) {
             L <- rlinegrid(angle, spacing, W)
             L <- L[W]
             do.call.matched("plot.psp",
                             resolve.defaults(list(x=L, add=TRUE),
                                              list(...)),
                             extrargs=c("lwd","lty","col"))
           }
         },
         polygonal = {
           p <- W$bdry
           # Determine whether user wants to fill the interior
           col.poly <- resolve.defaults(list(...), list(col=NA))$col
           den.poly <- resolve.defaults(list(...), list(density=NULL))$density
           no.fill  <- is.null(den.poly) &&
                       (is.null(col.poly) || is.na(col.poly))
           # Determine whether we need to triangulate the interior.
           # If it is required to fill the interior,
           # this can be done directly using polygon() provided
           # there are no holes. Otherwise we must triangulate the interior.
           if(no.fill)
             triangulate <- FALSE
           else {
             # Determine whether there are any holes
             holes <- (unlist(lapply(p, area.xypolygon)) < 0)
             triangulate <- any(holes)
           }

           if(!triangulate) {
             # No triangulation required;
             # simply plot the polygons
             for(i in seq(p))
               do.call.matched("polygon",
                               resolve.defaults(
                                                list(x=p[[i]]),
                                                list(...)),
                               extrargs="lwd")
           } else {
             if(!(spatstat.options("gpclib") && require(gpclib)))
               warning("Can't plot filled polygons: gpclib unavailable")
             else {
               # First fill the polygon's interior with colour
               # Use gpclib to triangulate
               txy <- triangulate(owin2gpc(W))
               ntri <- floor(nrow(txy)/3)
               # Fill triangles with colour
               for(i in seq(ntri)) {
                 ri <- 3 * (i - 1) + 1:3
                 do.call.matched("polygon",
                                 resolve.defaults(
                                                  list(x=txy[ri,]),
                                                  list(border=NA),
                                                  list(...)))
               }
             }
             # Now draw polygon boundaries
             for(i in seq(p))
               do.call.matched("polygon",
                               resolve.defaults(
                                                list(x=p[[i]]),
                                                list(density=0, col=NA),
                                                list(...)),
                               extrargs="lwd")
           }
           if(hatch) {
             L <- rlinegrid(angle, spacing, W)
             L <- L[W]
             do.call.matched("plot.psp",
                             resolve.defaults(list(x=L, add=TRUE),
                                              list(...)),
                             extrargs=c("lwd","lty","col"))
           }
         },
         mask = {
           # capture 'col' argument and ensure it's at least 2 values
           coldefault <- c(par("bg"), par("fg"))
           col <- resolve.defaults(
                                   list(...),
                                   spatstat.options("par.binary"),
                                   list(col=coldefault)
                                   )$col
           if(length(col) == 1) {
             col <- unique(c(par("bg"), col))
             if(length(col) == 1) 
               col <- c(par("fg"), col)
           }
           do.call.matched("image.default",
                           resolve.defaults(
                           list(x=W$xcol, y=W$yrow, z=t(W$m), add=TRUE),
                           list(col=col),       
                           list(...),
                           spatstat.options("par.binary")))
           if(hatch)
             warning("Hatching is not implemented for mask windows")
         },
         stop(paste("Don't know how to plot window of type", sQuote(W$type)))
         )

  invisible()
}





