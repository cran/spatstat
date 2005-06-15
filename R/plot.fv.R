#
#       plot.fv.R   (was: conspire.S)
#
#  $Revision: 1.13 $    $Date: 2005/06/14 04:27:42 $
#
#

conspire <- function(...) {
  .Deprecated("plot.fv")
  plot.fv(...)
}

plot.fv <- function(x, fmla, subset=NULL, lty=NULL, col=NULL, lwd=NULL,
                     xlim, ylim, xlab, ylab, ...) {

  verifyclass(x, "fv")
  indata <- as.data.frame(x)

  defaultplot <- missing(fmla)
  if(defaultplot)
    fmla <- attr(x, "fmla")

  # This *is* the last possible moment, so...
  fmla <- as.formula(fmla)

  # expand "."
  dotnames <- attr(x, "dotnames")
  if(is.null(dotnames)) {
    argu <- attr(x, "argu")
    allvars <- names(x)
    dotnames <- allvars[allvars != argu]
    dotnames <- rev(dotnames) # convention
  }
  u <- as.call(lapply(c("cbind", dotnames), as.name))
  fmla <- eval(substitute(substitute(fom, list(.=u)), list(fom=fmla)))

  # extract LHS and RHS
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

  # process lty, col, lwd arguments

  fixit <- function(a, n, a0) {
    if(is.null(a))
      return(a0)
    else if(length(a) == 1)
      return(rep(a, n))
    else if(length(a) != n)
      stop(paste("Length of", deparse(substitute(a)),
                 "does not match number of curves to be plotted"))
    else 
      return(a)
  }

  lty <- fixit(lty, nplots, 1:nplots)
  col <- fixit(col, nplots, 1:nplots)
  lwd <- fixit(lwd, nplots, rep(1, nplots))

  # plot lines
  
  for(i in 1:nplots)
    lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i], lwd=lwd[i])

  if(nplots == 1)
    return(invisible(NULL))
  else 
    return(data.frame(lty=lty, col=col, row.names=colnames(lhsdata)))
}


