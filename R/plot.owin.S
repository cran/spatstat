#
#	plot.owin.S
#
#	The 'plot' method for observation windows (class "owin")
#
#	$Revision: 1.9 $	$Date: 2002/07/18 10:18:45 $
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
      plot(x, y, type="n", ...,
           main=main, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i",
           asp=1.0)
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
      plot(x, y, xlim=xlim, ylim=ylim, type="n", ...,
           main=main, axes=FALSE, xlab="", ylab="", xaxs="i", yaxs="i")

    }
  }
# Draw window

  switch(W$type,
         rectangle = {
         },
         polygonal = {
           p <- W$bdry
           if(exists("is.R") && is.R()) {
             for(i in seq(p))
               polygon(p[[i]], ...)
           } else {
             for(i in seq(p))
               polygon(p[[i]], density=0, ...)
           }
         },
         mask = {
           # image(W$xcol, W$yrow, !t(W$m), add=TRUE, ...)
           arg1 <- list(W$xcol, W$yrow, !t(W$m), add=TRUE)
           arg2 <- spatstat.options("par.binary")[[1]]
           arg3 <- list(...)
           argue <- append(arg1, append(arg2, arg3))
           do.call("image", argue)
         },
         stop("Don't know how to plot image type", W$type)
         )

# Draw surrounding box
  wantbox <- (!missing(box) && box) || (missing(box) && W$type == "rectangle")
  if(wantbox)
    segments(x[c(1,2,2,1)],
             y[c(1,1,2,2)],
             x[c(2,2,1,1)],
             y[c(1,2,2,1)], ...)

  invisible()
}





