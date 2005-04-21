#
#	plot.ppp.S
#
#	$Revision: 1.17 $	$Date: 2005/04/20 20:37:37 $
#
#
#--------------------------------------------------------------------------

plot.ppp <-
  function(x, main, ..., chars, cols, use.marks=TRUE, add=FALSE, maxsize, markscale)
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

        
  if(!is.marked(x) || !use.marks) {
    do.call("points",
            resolve.defaults(list(x$x, x$y),
                             if(!missing(cols)) list(col=cols[1]) else NULL,
                             if(!missing(chars)) list(pch=chars[1]) else NULL,
                             list(...)))
    return(invisible())
  }

  # marked point pattern

  marks <- x$marks

  if(is.numeric(marks)) {
    # real-valued marks
    ok <- all(is.finite(marks))
    if(ok) {
      # examine spread of values
      mr <- range(marks)
      maxabs <- max(abs(mr))
      constant <- (diff(mr) < 4 * .Machine$double.eps)
      tiny <- (maxabs < 4 * .Machine$double.eps)
      if(tiny && missing(markscale)) {
        # data cannot be scaled successfully;
        # plot as points
        do.call("points",
                resolve.defaults(list(x$x, x$y),
                                 if(!missing(cols)) list(col=cols[1]) else NULL,
                                 list(...)))
        return(invisible())
      }
      # determine physical scale:
      #   plotted size = scal * (mark value) 
      if(!missing(maxsize) && !missing(markscale))
        stop("Only one of maxsize and markscale should be given")
      if(!missing(markscale))
        scal <- markscale
      else {
        # scale to be determined indirectly
        if(missing(maxsize)) {
          # guess appropriate max physical size of symbols
          maxsize <- 1.4/sqrt(pi * x$n/area.owin(x$window))
          maxsize <- min(maxsize, diameter(x$window) * 0.07)
        }
        # scale to [0,maxsize]
        scal <- maxsize/maxabs
      }
      # scale determined.
      # Apply the scaling
      ms <- marks * scal 

      # Finally, plot them..
      # plot positive values as circles
      neg <- (marks < 0)
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
  
  um <- if(is.factor(x$marks))
    levels(x$marks)
  else
    sort(unique(x$marks))

  ntypes <- length(um)
  
  if(missing(chars))
    chars <- seq(um)
  else if((nchars <- length(chars)) != ntypes) {
    if(nchars != 1)
      stop("length of \`chars\' is not equal to the number of types")
    else
      chars <- rep(chars, ntypes)
  }

  if(!missing(cols) && ((ncols <- length(cols)) != ntypes)) {
    if(ncols != 1)
      stop("length of \`cols\' is not equal to the number of types")
    else
      cols <- rep(cols, ntypes)
  }
    
  for(i in seq(um)) {
    relevant <- (x$marks == um[i])
    if(any(relevant))
      do.call("points",
              resolve.defaults(list(x$x[relevant], x$y[relevant]),
                               list(pch = chars[i]),
                               if(!missing(cols)) list(col=cols[i]) else NULL,
                               list(...)))
  }
  names(chars) <- um
  if(length(chars) < 20)
    return(chars)
  else
    return(invisible(chars))
}
