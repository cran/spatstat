#
# interactive plot for ppp objects using rpanel
#
#   $Revision: 1.8 $   $Date: 2010/03/08 08:23:04 $
#
#

# Effect:
# when the user types
#                 iplot(x)
# a pop-up panel displays a standard plot of x and
# buttons allowing control of the plot parameters.

# Coding:
# The panel 'p' contains the following internal variables
#      x          Original point pattern
#      xname      Name of x (for main title)
#      mtype      Type of marks of x
# The following variables in 'p' are controlled by panel buttons etc
#      split      Logical: whether to split multitype pattern
#      pointmap   Plot character, or "marks" indicating that marks are used
#      charsize   Character expansion factor cex
#      markscale  Mark scale factor markscale


iplot <- function(x, xname) {
  if(missing(xname))
    xname <- short.deparse(substitute(x))
  verifyclass(x, "ppp")
  if(markformat(x) == "dataframe")
	marks(x) <- marks(x)[,1]
  require(rpanel)
  mtype <- if(is.multitype(x)) "multitype" else if(is.marked(x)) "marked" else "unmarked"

  ##
  p <- rp.control(paste("iplot(", xname, ")", sep=""), 
                  x=x, xname=xname, mtype=mtype,
                  pointmap=if(is.marked(x)) "marks" else "o",
                  split=FALSE,
                  size=c(600, 400))

# Split panel into two halves  
# Left half of panel: display
# Right half of panel: controls
  rp.grid(p, "gdisplay", pos=list(row=0,column=0))
  rp.grid(p, "gcontrols", pos=list(row=0,column=1))

#----- Display side ------------

  # This line is to placate the package checker
  mytkr <- NULL
  
  rp.tkrplot(p, mytkr, do.iplot, pos=list(row=0,column=0,grid="gdisplay"))

  redraw <- function(panel) {
    rp.tkrreplot(p, mytkr)
    panel
  }
  
#----- Control side ------------
  nextrow <- 0
  pozzie <- function(n=nextrow) list(row=n,column=0,grid="gcontrols")
# main title
  rp.textentry(p, xname, action=redraw, title="Plot title",
               pos=pozzie(0))
  nextrow <- 1

# split ?
  if(mtype == "multitype") {
    rp.checkbox(p, split, initval=FALSE, 
                title="Split according to marks", action=redraw,
                pos=pozzie(1))
    nextrow <- 2
  }

# plot character or mark style
  ptvalues <- c("o", "bullet", "plus")
  ptlabels <- c("open circles", "filled circles", "crosshairs")
  if(is.marked(x)) {
    ptvalues <- c("marks", ptvalues)
    ptlabels <- if(mtype == "multitype")
      c("Symbols depending on mark", ptlabels)
    else c("Circles proportional to mark", ptlabels)
  }
  pointmap <- ptvalues[1]
  rp.radiogroup(p, pointmap, values=ptvalues, labels=ptlabels,
   			  title="how to plot points", action=redraw,
                pos=pozzie(nextrow))
  nextrow <- nextrow+1

# plot character size
  charsize <- 1
  rp.slider(p, charsize, 0, 5, action=redraw, 
            title="symbol expansion factor (cex)", initval=1, showvalue=TRUE,
            pos=pozzie(nextrow))
  nextrow <- nextrow+1
  
# mark scale
  if(mtype == "marked") {
    marx <- x$marks
    marx <- marx[is.finite(marx)]
    scal <- mark.scale.default(marx, x$window)
    markscale <- scal
    rp.slider(p, markscale, from=scal/10, to = 10*scal, action=redraw,
              initval=scal,
              title="mark scale factor (markscale)", showvalue=TRUE,
              pos=pozzie(nextrow))
    nextrow <- nextrow+1
  }

# button to print a summary at console
  rp.button(p, title="Print summary information",
            pos=pozzie(nextrow),
            action=function(panel) { print(summary(panel$x)); panel} )
  nextrow <- nextrow+1
# quit button 
  rp.button(p, title="Quit", quitbutton=TRUE, pos=pozzie(nextrow),
            action= function(panel) { panel })
#  
  invisible(NULL)
}

# function that updates the plot when the control panel is operated

do.iplot <- function(panel) { 
  use.marks <- TRUE
  pch <- 16
  switch(panel$pointmap,
         marks={
           use.marks <- TRUE
           pch <- NULL
         }, 
         o = {
           use.marks <- FALSE
           pch <- 1
         }, 
         bullet = {
           use.marks <- FALSE
           pch <- 16
         },
         plus = {
           use.marks <- FALSE
           pch <- 3
         })
  y <- panel$x
  if(panel$mtype == "multitype" && panel$split)
    y <- split(y, un=(panel$pointmap != "marks"))
  if(panel$mtype == "marked" && panel$pointmap == "marks") 
    plot(y, main=panel$xname, use.marks=use.marks, markscale=panel$markscale)
  else
    plot(y, main=panel$xname, use.marks=use.marks, 
         pch=pch, cex=panel$charsize)      
  panel
}

