#
#	plot.ppp.S
#
#	$Revision: 1.12 $	$Date: 2002/05/13 12:41:10 $
#
#
#--------------------------------------------------------------------------

plot.ppp <-
  function(x, main, ..., chars, cols, use.marks=TRUE, add=FALSE, maxsize)
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

        
  if(is.null(x$marks) || !use.marks) {
    points(x$x, x$y, ...)
    return(invisible())
  }

  # marked point pattern

  marks <- x$marks

  if(is.numeric(marks)) {
    # real-valued marks
    ok <- !(any(is.na(marks)) || any(is.infinite(marks)))
    if(ok) {
      # rescale marks to unit interval
      mr <- range(c(0,marks))
      if(diff(mr) > 4 * .Machine$double.eps)
        ms <- (marks - mr[1])/diff(mr)
      else
        ms <- rep(1, length(marks))
      # determine roughly appropriate max radius on plot
      if(missing(maxsize)) {
        maxsize <- 1.4/sqrt(pi * x$n/area.owin(x$window))
        maxsize <- min(maxsize, diameter(x$window) * 0.07)
      }
      # fire
      symbols(x$x, x$y, circles=ms * maxsize, inches=FALSE, add=TRUE, ...)
      # compute plottable scale bar for marks
      mp.value <- pretty(mr)
      mp.plotted <- mp.value * maxsize/mr[2]
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

  if(missing(chars))
    chars <- seq(um)

  # generic argument 'col' conflicts with 'cols'
  col.given <- ("col" %in% names(list(...)))
  cols.given <- !missing(cols)
  if(cols.given && col.given)
    stop("Only one of the arguments \"col\" and \"cols\" should be given")
  
  for(i in seq(um)) {
    relevant <- (x$marks == um[i])
    if(any(relevant)) {
      if(cols.given)
        points(x$x[relevant], x$y[relevant], pch = chars[i], col=cols[i], ...)
      else
        points(x$x[relevant], x$y[relevant], pch = chars[i], ...)
    }
  }
  names(chars) <- um
  if(length(chars) < 20)
    return(chars)
  else
    return(invisible(chars))
}
