#
#	plot.ppp.S
#
#	$Revision: 1.23 $	$Date: 2006/10/24 03:00:39 $
#
#
#--------------------------------------------------------------------------

plot.ppp <-
  function(x, main, ..., chars, cols, use.marks=TRUE, add=FALSE,
           maxsize=NULL, markscale=NULL)
{
#
# Function plot.ppp.
# A plot() method for the class 'ppp'
#
  if(missing(main))
    main <- deparse(substitute(x))

  x <- as.ppp(x)
  
  if(!add)
    plot.owin(x$window, ..., main=main)

  if(x$n == 0)
    return(invisible())

  # handle plot parameters
  explicit <- list()
  if(!missing(cols))
    explicit <- append(explicit, list(cols=cols))
  if(!missing(chars))
    explicit <- append(explicit, list(chars=chars))
    
  defaults <- spatstat.options("par.points")

  smartpoints <- function(..., index=1, col=NULL, pch=NULL, cols, chars) {
    if(missing(col) && !missing(cols))
      col <- cols[index]
    if(missing(pch) && !missing(chars))
      pch <- chars[index]
    do.call("points",
            resolve.defaults(list(...),
                             if(!is.null(col)) list(col=col) else NULL,
                             if(!is.null(pch)) list(pch=pch) else NULL))
  }

  marked <- is.marked(x, dfok=TRUE)
  if(!marked || !use.marks) {
    do.call("smartpoints",
            resolve.defaults(list(x$x, x$y),
                             explicit,
                             list(...),
                             spatstat.options("par.points")))
    return(invisible())
  }

  # marked point pattern

  marx <- marks(x, dfok=TRUE)

  # if it's a data frame, take the first column
  if(is.data.frame(marx))
    marx <- marx[,1]

  if(is.numeric(marx)) {
    # real-valued marks
    ok <- all(is.finite(marx))
    if(ok) {
      # establish values of markscale, maxsize
      if(!is.null(maxsize) && !is.null(markscale))
        stop("Only one of maxsize and markscale should be given")
      if(is.null(maxsize) && is.null(markscale)) {
        # if BOTH are absent, enforce the spatstat defaults
        # (which could also be null)
        pop <- spatstat.options("par.points")
        markscale <- pop$markscale
        maxsize   <- pop$maxsize
      }
      # examine spread of values
      mr <- range(marx)
      maxabs <- max(abs(mr))
      constant <- (diff(mr) < 4 * .Machine$double.eps)
      tiny <- (maxabs < 4 * .Machine$double.eps)
      if(tiny && is.null(markscale)) {
        # data cannot be scaled successfully;
        # plot as points
        do.call("smartpoints",
                resolve.defaults(list(x$x, x$y),
                                 explicit,
                                 list(...),
                                 spatstat.options("par.points")))
        return(invisible())
      }
      # determine physical scale:
      #   plotted size = scal * (mark value) 
      if(!is.null(markscale))
        scal <- markscale
      else {
        # scale to be determined indirectly
        if(is.null(maxsize)) {
          # guess appropriate max physical size of symbols
          maxsize <- 1.4/sqrt(pi * x$n/area.owin(x$window))
          maxsize <- min(maxsize, diameter(x$window) * 0.07)
        }
        # scale to [0,maxsize]
        scal <- maxsize/maxabs
      }
      # scale determined.
      # Apply the scaling
      ms <- marx * scal 

      # Finally, plot them..
      # plot positive values as circles
      neg <- (marx < 0)
      if(any(!neg)) 
        do.call("symbols",
                resolve.defaults(
                                 list(x$x[!neg], x$y[!neg]),
                                 list(circles = ms[!neg]),
                                 list(inches = FALSE, add = TRUE),
                                 if(!missing(cols)) list(fg=cols[1]) else NULL,
                                 list(...)))
      # plot negative values as squares
      if(any(neg))
        do.call("symbols",
                resolve.defaults(
                                 list(x$x[neg], x$y[neg]),
                                 list(squares = - ms[neg]),
                                 list(inches = FALSE, add = TRUE),
                                 if(!missing(cols)) list(fg=cols[1]) else NULL,
                                 list(...)))
      # return a plottable scale bar
      mp.value <- if(constant) mr[1] else pretty(mr)
      mp.plotted <- mp.value * scal 
      names(mp.plotted) <- paste(mp.value)
      return(mp.plotted)
    } else {
      warning("Some marks are NA or Inf; treating marks as non-numeric")
    }
  }
  
  um <- if(is.factor(marx))
    levels(marx)
  else
    sort(unique(marx))

  ntypes <- length(um)
  
  if(missing(chars))
    chars <- seq(um)
  else if((nchars <- length(chars)) != ntypes) {
    if(nchars != 1)
      stop(paste("length of",
                 sQuote("chars"),
                 "is not equal to the number of types"))
    else
      chars <- rep(chars, ntypes)
  }

  if(!missing(cols) && ((ncols <- length(cols)) != ntypes)) {
    if(ncols != 1)
      stop(paste("length of",
                 sQuote("cols"),
                 "is not equal to the number of types"))
    else
      cols <- rep(cols, ntypes)
  }
    
  for(i in seq(um)) {
    relevant <- (marx == um[i])
    if(any(relevant))
      do.call("smartpoints",
              resolve.defaults(list(x$x[relevant], x$y[relevant]),
                               list(pch = chars[i]),
                               explicit,
                               list(index=i),
                               list(...),
                               spatstat.options("par.points")))
  }
  names(chars) <- um
  if(length(chars) < 20)
    return(chars)
  else
    return(invisible(chars))
}
