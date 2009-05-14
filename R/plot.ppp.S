#
#	plot.ppp.S
#
#	$Revision: 1.37 $	$Date: 2009/05/12 02:12:05 $
#
#
#--------------------------------------------------------------------------

plot.ppp <-
  function(x, main, ..., chars=NULL, cols=NULL, use.marks=TRUE, add=FALSE,
           maxsize=NULL, markscale=NULL)
{
#
# Function plot.ppp.
# A plot() method for the class 'ppp'
#
  if(missing(main))
    main <- deparse(substitute(x))

# First handle `rejected' points
  sick <- inherits(x, "ppp") && !is.null(rejects <- attr(x, "rejects"))
  if(sick) {
    # get any parameters
    par.direct <- list(main=main, use.marks=use.marks,
                   maxsize=maxsize, markscale=markscale)
    par.rejects.default <- list(pch="+")
    par.rejects <- resolve.defaults(list(...),
                                    list(par.rejects=par.rejects.default))$par.rejects
    par.rejects <- resolve.defaults(par.rejects, par.rejects.default)
    par.all <- resolve.defaults(par.rejects, par.direct)
    rw <- resolve.defaults(list(...), list(rejectwindow=NULL))$rejectwindow
    # determine window for rejects
    rwin <-
      if(is.null(rw))
        rejects$window
      else if(is.logical(rw) && rw)
        rejects$window
      else if(inherits(rw, "owin"))
        rw
      else if(is.character(rw)) {
        switch(rw,
               box={bounding.box(rejects, x)},
               ripras={ripras(c(rejects$x, x$x), c(rejects$y, x$y))},
               stop(paste("Unrecognised option: rejectwindow=", rw)))
      } else stop("Unrecognised format for rejectwindow")
    if(is.null(rwin))
      stop("Selected window for rejects pattern is NULL")
    # Create suitable space
    plot(rejects$window, add=add, type="n", main="")
    if(!add)
      title(main=main)
    # plot window if commanded
    if(!is.null(rw)) {
      rwinpardefault <- list(lty=2,lwd=1,border=1)
      rwinpars <-
        resolve.defaults(par.rejects, rwinpardefault)[names(rwinpardefault)]
      do.call("plot.owin", append(list(rwin, add=TRUE), rwinpars))
    }
    # plot points
    do.call("plot.ppp", append(list(rejects, add=TRUE), par.all))
    plot.owin(x$window, add=TRUE, ...)
    warning(paste(rejects$n, "illegal points also plotted"))
    # the rest is added
    add <- TRUE
  }

# Now convert to bona fide point pattern
  x <- as.ppp(x)
  
  if(!add)
    plot.owin(x$window, ..., main=main)
    
  if(x$n == 0)
    return(invisible())

  # handle plot parameters
  explicit <- list()
  if(!is.null(cols))
    explicit <- append(explicit, list(cols=cols))
  if(!is.null(chars))
    explicit <- append(explicit, list(chars=chars))
    
  defaults <- spatstat.options("par.points")

  smartpoints <- function(xx, yy, ...,
                          index=1, col=NULL, pch=NULL, cols=NULL, chars=NULL) {
    if(is.null(col) && !is.null(cols))
      col <- cols[index]
    if(is.null(pch) && !is.null(chars))
      pch <- chars[index]
    do.call.matched("points",
            resolve.defaults(list(x=list(x=xx, y=yy), ...),
                             if(!is.null(col)) list(col=col) else NULL,
                             if(!is.null(pch)) list(pch=pch) else NULL),
                    extrargs=c("col", "pch", "type", "bg", "cex", "lwd"))
  }

  marked <- is.marked(x, dfok=TRUE)
  if(!marked || !use.marks) {
    do.call("smartpoints",
            resolve.defaults(list(xx=x$x, yy=x$y),
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

  # check there are some valid marks!
  ok <- !is.na(marx)
  if(all(!ok)) {
    warning("All mark values are NA; plotting locations only.")
    do.call("smartpoints",
            resolve.defaults(list(xx=x$x, yy=x$y),
                             explicit,
                             list(...),
                             spatstat.options("par.points")))
    return(invisible())
  }

  # otherwise ignore invalid marks
  if(!all(ok)) {
    warning(paste("Some marks are NA;",
                    "corresponding points are omitted."))
    x <- x[ok]
    marx <- marx[ok]
  }

  ################  real-valued marks ############################

  if(is.numeric(marx)) {

    ok <- is.finite(marx)

    if(!all(ok)) {
      warning(paste("Some marks are infinite",
                    "corresponding points are omitted."))
      x <- x[ok]
      marx <- marx[ok]
    }

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
        bb <- as.rectangle(x$window)
        maxsize <- 1.4/sqrt(pi * x$n/area.owin(bb))
        maxsize <- min(maxsize, diameter(bb) * 0.07)
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
                               if(!is.null(cols)) list(fg=cols[1]) else NULL,
                               list(...)))
    # plot negative values as squares
    if(any(neg))
      do.call("symbols",
              resolve.defaults(
                               list(x$x[neg], x$y[neg]),
                               list(squares = - ms[neg]),
                               list(inches = FALSE, add = TRUE),
                               if(!is.null(cols)) list(fg=cols[1]) else NULL,
                               list(...)))
    # return a plottable scale bar
    mp.value <- if(constant) mr[1] else pretty(mr)
    mp.plotted <- mp.value * scal 
    names(mp.plotted) <- paste(mp.value)
    return(mp.plotted)
  }

  ##################### non-numeric marks ###############################

  um <- if(is.factor(marx))
    levels(marx)
  else
    sort(unique(marx))

  ntypes <- length(um)
  
  if(is.null(chars)) {
    if(ntypes <= 25) {
      # numerical 'pch' 
      chars <- 1:ntypes
    } else if(ntypes <= 52) {
      # letters
      ltr <- c(letters, LETTERS)
      chars <- ltr[1:ntypes]
    } else {
      # wrapped sequence of letters
      warning("There are too many types to display every type as a different character")
      chars <- ltr[1 + (0:(ntypes - 1) %% 52)]
    }
  }
  else if((nchars <- length(chars)) != ntypes) {
    if(nchars != 1)
      stop(paste("length of",
                 sQuote("chars"),
                 "is not equal to the number of types"))
    else
      chars <- rep(chars, ntypes)
  }

  if(!is.null(cols) && ((ncols <- length(cols)) != ntypes)) {
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

