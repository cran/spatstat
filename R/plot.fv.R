#
#       plot.fv.R   (was: conspire.S)
#
#  $Revision: 1.8 $    $Date: 2005/02/08 01:52:32 $
#
#

conspire <- function(...) {
  .Deprecated("plot.fv")
  plot.fv(...)
}

plot.fv <- function(x, fmla, subset=NULL, lty=NULL, col=NULL,
                     xlim, ylim, xlab, ylab, ...) {

  verifyclass(x, "fv")
  indata <- as.data.frame(x)

  defaultplot <- missing(fmla)
  if(defaultplot)
    fmla <- attr(x, "fmla")

  # This *is* the last possible moment, so...
  fmla <- as.formula(fmla)

  lhs <- fmla[[2]]
  rhs <- fmla[[3]]

  # evaluate expression a in data frame b
  evaluate <- function(a,b) {
    if(exists("is.R") && is.R())
      eval(a, envir=b)
    else
      eval(a, local=b)
  }
  
  lhsdata <- evaluate(lhs, indata)
  rhsdata <- evaluate(rhs, indata)
  
  if(is.vector(lhsdata))
    lhsdata <- matrix(lhsdata, ncol=1)

  if(!is.vector(rhsdata))
    stop("rhs of formula seems not to be a vector")

  # restrict data to subset if desired
  if(!is.null(subset)) {
    keep <- if(is.character(subset))
		evaluate(parse(text=subset), indata)
            else
                evaluate(subset, indata)
    lhsdata <- lhsdata[keep, , drop=FALSE]
    rhsdata <- rhsdata[keep]
  } 

  # determine x and y limits and clip data to these limits
  if(!missing(xlim)) {
    ok <- (xlim[1] <= rhsdata & rhsdata <= xlim[2])
    rhsdata <- rhsdata[ok]
    lhsdata <- lhsdata[ok, , drop=FALSE]
  } else {
    # if we're using the default argument, use its recommended range
    if(rhs == attr(x, "argu")) {
      xlim <- attr(x,"alim")
      ok <- is.finite(rhsdata) & rhsdata >= xlim[1] & rhsdata <= xlim[2]
      rhsdata <- rhsdata[ok]
      lhsdata <- lhsdata[ok, , drop=FALSE]
    } else { # actual range of values to be plotted
      xlim <- range(rhsdata[is.finite(rhsdata)],na.rm=TRUE)
      rok <- is.finite(rhsdata) & rhsdata >= xlim[1] & rhsdata <= xlim[2]
      lok <- apply(is.finite(lhsdata), 1, any)
      ok <- lok & rok
      rhsdata <- rhsdata[ok]
      lhsdata <- lhsdata[ok, , drop=FALSE]
      xlim <- range(rhsdata)
    }
  }
  
  if(missing(ylim))
    ylim <- range(lhsdata,na.rm=TRUE)

  # work out how to label the plot
  if(missing(xlab))
    xlab <- as.character(fmla)[3]

  if(missing(ylab)) {
    yl <- attr(x, "ylab")
    if(!is.null(yl) && defaultplot)
      ylab <- yl
    else {
      yname <- paste(lhs)
      if(length(yname) > 1 && yname[[1]] == "cbind")
        ylab <- paste(yname[-1], collapse=" , ")
      else
        ylab <- as.character(fmla)[2]
    }
  }

  # check for argument "add"=TRUE
  dotargs <- list(...)
  v <- match("add", names(dotargs))
  if(is.na(v) || !(addit <- as.logical(dotargs[[v]])))
    plot(xlim, ylim, type="n", xlab=xlab, ylab=ylab, ...)

  nplots <- ncol(lhsdata)

  if(is.null(lty))
    lty <- 1:nplots
  else if(length(lty) == 1)
    lty <- rep(lty, nplots)
  else if(length(lty) != nplots)
    stop("Length of \`lty\' does not match number of curves to be plotted")
  
  if(is.null(col))
    col <- 1:nplots
  else if(length(col) == 1)
    col <- rep(col, nplots)
  else if(length(col) != nplots)
    stop("Length of \`col\' does not match number of curves to be plotted")
  
  for(i in 1:nplots)
    lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i])

  if(nplots == 1)
    return(invisible(NULL))
  else 
    return(data.frame(lty=lty, col=col, row.names=colnames(lhsdata)))
}


