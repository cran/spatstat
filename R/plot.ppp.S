#
#	plot.ppp.S
#
#	$Revision: 1.13 $	$Date: 2003/03/11 02:55:10 $
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
      # guess appropriate max physical size of symbols
      if (missing(maxsize)) {
        maxsize <- 1.4/sqrt(pi * x$n/area.owin(x$window))
        maxsize <- min(maxsize, diameter(x$window) * 0.07)
      }
      # find range of values
      mr <- range(c(0,marks))
      maxabs <- max(abs(mr))
      # determine physical scale and apply it
      if (diff(mr) < 4 * .Machine$double.eps
          || maxabs < 4 * .Machine$double.eps)
        {
          # constant values - plot at half maxsize
          ms <- rep(0.5 * maxsize, length(marks))
          mp.value <- mr[1]
          mp.plotted <- 0.5 * maxsize
        } else {
          # scale to [0,maxsize]
          scal <- maxsize/maxabs
          ms <- marks * scal
          mp.value <- pretty(mr)
          mp.plotted <- mp.value * scal
        }
      # plot positive values as circles
      neg <- (marks < 0)
      if(any(!neg))
        symbols(x$x[!neg], x$y[!neg],
                circles = ms[!neg],
                inches = FALSE, add = TRUE, ...)
      # plot negative values as squares
      if(any(neg))
        symbols(x$x[neg], x$y[neg],
                squares = - ms[neg],
                inches = FALSE, add = TRUE, ...)

      # return a plottable scale bar
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
