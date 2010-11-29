#
#       plot.fv.R   (was: conspire.S)
#
#  $Revision: 1.50 $    $Date: 2010/11/22 04:40:28 $
#
#

conspire <- function(...) {
  .Deprecated("plot.fv", package="spatstat")
  plot.fv(...)
}

plot.fv <- function(x, fmla, ..., subset=NULL, lty=NULL, col=NULL, lwd=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                    ylim.covers=NULL, legend=TRUE, legendpos="topleft",
                    legendmath=FALSE, shade=NULL, shadecol="grey") {

  xname <-
    if(is.language(substitute(x))) deparse(substitute(x)) else ""
  verifyclass(x, "fv")

  indata <- as.data.frame(x)

  defaultplot <- missing(fmla) || is.null(fmla)
  if(defaultplot) 
    fmla <- attr(x, "fmla")

  # This *is* the last possible moment, so...
  fmla <- as.formula(fmla)

  # Extract left hand side as given
  lhs.original <- fmla[[2]]
  
  # expand "."
  dotnames <- fvnames(x, ".")
  if(is.null(dotnames)) {
    argu <- fvnames(x, ".x")
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
#    if(exists("is.R") && is.R())
      eval(a, envir=b)
#    else
#      eval(a, local=b)
  }
  
  lhsdata <- evaluate(lhs, indata)
  rhsdata <- evaluate(rhs, indata)
  
  if(is.vector(lhsdata)) {
    lhsdata <- matrix(lhsdata, ncol=1)
    colnames(lhsdata) <- paste(lhs, collapse="")
  }

  if(is.matrix(rhsdata))
    stop("rhs of formula should yield a vector")
  rhsdata <- as.numeric(rhsdata)

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
  if(!is.null(xlim)) {
    ok <- (xlim[1] <= rhsdata & rhsdata <= xlim[2])
    rhsdata <- rhsdata[ok]
    lhsdata <- lhsdata[ok, , drop=FALSE]
  } else {
    # if we're using the default argument, use its recommended range
    if(rhs == fvnames(x, ".x")) {
      xlim <- attr(x, "alim")
      rok <- is.finite(rhsdata) & rhsdata >= xlim[1] & rhsdata <= xlim[2]
      lok <- apply(is.finite(lhsdata), 1, any)
      ok <- lok & rok
      rhsdata <- rhsdata[ok]
      lhsdata <- lhsdata[ok, , drop=FALSE]
    } else { # actual range of values to be plotted
      rok <- is.finite(rhsdata) 
      lok <- apply(is.finite(lhsdata), 1, any)
      ok <- lok & rok
      rhsdata <- rhsdata[ok]
      lhsdata <- lhsdata[ok, , drop=FALSE]
      xlim <- range(rhsdata)
    }
  }
  
  if(is.null(ylim))
    ylim <- range(lhsdata[is.finite(lhsdata)],na.rm=TRUE)
  if(!is.null(ylim.covers))
    ylim <- range(ylim, ylim.covers)

  # work out how to label the plot
  
  # extract plot labels 
  labl <- attr(x, "labl")
  # expand plot labels
  if(!is.null(fname <- attr(x, "fname")))
    labl <- sprintf(labl, fname)
  # construct mapping from identifiers to labels
  map <- as.list(labl)
  magic <- function(x) {
    subx <- paste("substitute(", x, ", NULL)")
    eval(parse(text=subx))
  }
  map <- lapply(map, magic)
  names(map) <- colnames(x)
  # also map "." to name of target function
  if(!is.null(ye <- attr(x, "yexp")))
        map <- append(map, list("."=ye))
  # alternative version of map (vector of expressions)
  mapvec <- sapply(as.list(labl), function(x) { parse(text=x) })
  names(mapvec) <- colnames(x)

  # ......... label for x axis ..................

  if(is.null(xlab)) {
    argname <- fvnames(x, ".x")
    if(as.character(fmla)[3] == argname) {
      # the x axis is the default function argument.
      # Add name of unit of length
      ax <- summary(unitname(x))$axis
      xlab <- paste(argname, ax)
    } else {
      # map ident to label
      xlab <- eval(substitute(substitute(rh, mp), list(rh=rhs, mp=map)))
    }
  }
  if(is.language(xlab) && !is.expression(xlab))
    xlab <- as.expression(xlab)
      
  # ......... label for y axis ...................
  
  if(is.null(ylab)) {
    yl <- attr(x, "yexp")
    if(defaultplot && !is.null(yl)) {
      ylab <- yl
    } else {
      # make y axis label from original LHS
      leftside <- lhs.original
      if(ncol(lhsdata) > 1) {
        # replace 'cbind(....)' by '.' for labelling purposes only
        leftside <- paste(as.expression(leftside))
        cb <- paste("cbind(",
                    paste(colnames(lhsdata), collapse=", "),
                    ")", sep="")
        leftside <- gsub(cb, ".", leftside, fixed=TRUE)
        # convert back to language
        leftside <- as.formula(paste(leftside, "~1"))[[2]]
      }
      # replace short identifiers by plot labels
      ylab <- eval(substitute(substitute(le, mp),
                                list(le=leftside, mp=map)))
    }
  }
  if(is.language(ylab) && !is.expression(ylab))
    ylab <- as.expression(ylab)

  # check for argument "add" (defaults to FALSE)
  dotargs <- list(...)
  v <- match("add", names(dotargs))
  addit <- if(is.na(v)) FALSE else as.logical(dotargs[[v]])
  
  # if add=FALSE, create new plot
  if(!addit)
    do.call("plot.default",
            resolve.defaults(list(xlim, ylim, type="n"),
                             list(xlab=xlab, ylab=ylab),
                             list(...),
                             list(main=xname)))

  nplots <- ncol(lhsdata)

  # process lty, col, lwd arguments

  fixit <- function(a, n, a0, a00) {
    if(is.null(a))
      a <- if(!is.null(a0)) a0 else a00
    if(length(a) == 1)
      return(rep(a, n))
    else if(length(a) != n)
      stop(paste("Length of", deparse(substitute(a)),
                 "does not match number of curves to be plotted"))
    else 
      return(a)
  }

  opt0 <- spatstat.options("par.fv")
  
  lty <- fixit(lty, nplots, opt0$lty, 1:nplots)
  col <- fixit(col, nplots, opt0$col, 1:nplots)
  lwd <- fixit(lwd, nplots, opt0$lwd, 1)

  allind <- 1:nplots

  if(!is.null(shade)) {
    # shade region between critical boundaries
    # select columns by name or number
    names(allind) <- colnames(lhsdata)
    shind <- try(allind[shade])
    if(inherits(shind, "try-error")) 
      stop("The argument shade should be a valid subset index for columns of x")
    if(any(nbg <- is.na(shind))) {
      # columns not included in formula; get them
      morelhs <- try(as.matrix(indata[ok, shade[nbg], drop=FALSE]))
      if(inherits(morelhs, "try-error")) 
        stop("The argument shade should be a valid subset index for columns of x")
      nmore <- ncol(morelhs)
      lhsdata <- cbind(lhsdata, morelhs)
      shind[nbg] <- nplots + seq(nmore)
      nplots <- nplots + nmore
      lty <- c(lty, rep(lty[1], nmore))
      col <- c(col, rep(col[1], nmore))
      lwd <- c(lwd, rep(lwd[1], nmore))
    }
    # extract relevant columns
    shdata <- lhsdata[, shind]
    if(!is.matrix(shdata) || ncol(shdata) != 2) 
      stop("The argument shade should select two columns of x")
    # plot grey polygon between these limits
    shadeOK <- is.finite(rhsdata) & apply(is.finite(shdata), 1, all)
    polygon(c(rhsdata[shadeOK], rev(rhsdata[shadeOK])),
            c(shdata[shadeOK,1],  rev(shdata[shadeOK,2])),
            border=shadecol, col=shadecol)
    # overwrite graphical parameters
    lty[shind] <- 1
    # try to preserve the same type of colour specification
    if(is.character(col) && is.character(shadecol)) {
      # character representations 
      col[shind] <- shadecol
    } else if(is.numeric(col) && !is.na(sc <- paletteindex(shadecol))) {
      # indices in colour palette
      col[shind] <- sc
    } else {
      # convert colours to hexadecimal and edit relevant values
      col <- col2hex(col)
      col[shind] <- col2hex(shadecol)
    }
    # remove these columns from further plotting
    allind <- allind[-shind]
    # 
  }
  
  # plot lines

  for(i in allind)
    lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i], lwd=lwd[i])

  if(nplots == 1)
    return(invisible(NULL))
  else {
    key <- colnames(lhsdata)
    mat <- match(key, names(x))
    desc <- attr(x, "desc")[mat]
    labl <- labl[mat]
    ylab <- attr(x, "ylab")
    if(!is.null(ylab)) {
      if(is.language(ylab))
        ylab <- deparse(ylab)
      desc <- sprintf(desc, ylab)
    }
    if(!is.null(legend) && legend) {
      # do legend
      legtxt <- key
      if(legendmath) {
        legtxt <- labl
        # try to convert to expressions
        fancy <- try(parse(text=labl))
        if(!inherits(fancy, "try-error"))
          legtxt <- fancy
      }
      legend(legendpos, inset=0.05, lty=lty, col=col, legend=legtxt)
    }
    df <- data.frame(lty=lty, col=col, key=key, label=labl,
                      meaning=desc, row.names=key)
    return(df)
  }
}

