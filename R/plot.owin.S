#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.15 $	$Date: 2006/06/05 09:32:53 $
#
#
#

plot.owin <- function(x, main, add=FALSE, ..., box=TRUE, edge=0.04)
{
#
# Function plot.owin.  A method for plot.
#
# argument must be called 'x' for compatibility with plot()
  W <- x
  if(missing(main))
    main <- deparse(substitute(x))
# no, this cannot be inserted in the argument list!!!
  verifyclass(W, "owin")

#########        
  x <- W$xrange
  y <- W$yrange

####################################################  
  if(!add) {
# create new plotting region    
# check whether argument 'asp' is recognised
    plargs <- names(formals(plot.default))
    if(!is.null(plargs) && "asp" %in% plargs) {
# set up plot with equal scales
      do.call.matched("plot.default",
                      resolve.defaults(list(x=x, y=y, type="n"),
                                       list(...),
                                       list(main=main, axes=FALSE,
                                            xlab="", ylab="",
                                            xaxs="i", yaxs="i",
                                            asp=1.0)))
    } else {
# D.I.Y.            
# see 'help(par)' under 'xaxs': default is xaxs="r";
# we seize full control by setting xaxs="i",
# but mimic the extra 4% space allocated when xaxs="r"
        
      blowup <- function(v, s) { mean(v) + s * (v - mean(v)) }
      xlim <- blowup(x, 1+edge)
      ylim <- blowup(y, 1+edge)

# Constrain X scale = Y scale        
      pin <- par("pin")  #physical size of plot region
      xscale <- pin[1]/diff(xlim)
      yscale <- pin[2]/diff(ylim)
      sc <- min(xscale,yscale)
      if(xscale > sc) xlim <- blowup(xlim, xscale/sc)
      if(yscale > sc) ylim <- blowup(ylim, yscale/sc) 

# Commit scales
      do.call.matched("plot.default",
                      resolve.defaults(list(x=x, y=y, type="n"),
                                       list(...),
                                       list(xlim=xlim, ylim=ylim, main=main,
                                            axes=FALSE, xlab="", ylab="",
                                            xaxs="i", yaxs="i")))
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
                                     list(x0=x[c(1,2,2,1)],
                                          y0=y[c(1,1,2,2)],
                                          x1=x[c(2,2,1,1)],
                                          y1=y[c(1,2,2,1)]),
                                     list(...)))
  invisible()
}





