#
#       plot.fv.R   (was: conspire.S)
#
#  $Revision: 1.74 $    $Date: 2011/07/06 03:38:36 $
#
#

conspire <- function(...) {
  .Deprecated("plot.fv", package="spatstat")
  plot.fv(...)
}

plot.fv <- function(x, fmla, ..., subset=NULL, lty=NULL, col=NULL, lwd=NULL,
                    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
                    ylim.covers=NULL, legend=!add, legendpos="topleft",
                    legendmath=TRUE, legendargs=list(),
                    shade=NULL, shadecol="grey", add=FALSE) {

  xname <-
    if(is.language(substitute(x))) deparse(substitute(x)) else ""

  verifyclass(x, "fv")
  env.user <- parent.frame()

  indata <- as.data.frame(x)

  # ---------------- determine plot formula ----------------
  
  defaultplot <- missing(fmla) || is.null(fmla)
  if(defaultplot) 
    fmla <- attr(x, "fmla")

  # This *is* the last possible moment, so...
  fmla <- as.formula(fmla, env=env.user)

  # validate the variable names
  vars <- variablesinformula(fmla)
  reserved <- c(".", ".x", ".y")
  external <- !(vars %in% c(colnames(x), reserved))
  if(any(external)) {
    sought <- vars[external]
    found <- unlist(lapply(sought, exists, envir=env.user))
    if(any(!found)) {
      nnot <- sum(!found)
      stop(paste(ngettext(nnot, "Variable", "Variables"),
                 commasep(sQuote(sought[!found])),
                 ngettext(nnot, "was", "were"),
                 "not found"))
    } else {
      # validate the found variables
      externvars <- lapply(sought, get, envir=env.user)
      ok <- unlist(lapply(externvars,
                          function(z, n) { is.numeric(z) &&
                                           length(z) %in% c(1,n) },
                          n=nrow(x)))
      if(!all(ok)) {
        nnot <- sum(!ok)
        stop(paste(ngettext(nnot, "Variable", "Variables"),
                   commasep(sQuote(sought[!ok])),
                   ngettext(nnot, "is", "are"),
                   "not of the right format"))
      }
    }
  }
  
  # Extract left hand side as given
  lhs.original <- fmla[[2]]
  
  # expand "."
  dotnames <- fvnames(x, ".")
  u <- as.call(lapply(c("cbind", dotnames), as.name))
  ux <- as.name(fvnames(x, ".x"))
  uy <- as.name(fvnames(x, ".y"))
  fmla <- eval(substitute(substitute(fom, list(.=u, .x=ux, .y=uy)),
                          list(fom=fmla)))

  # ------------------- extract data for plot ---------------------
  
  # extract LHS and RHS of formula
  lhs <- fmla[[2]]
  rhs <- fmla[[3]]

  # extract data 
  lhsdata <- eval(lhs, envir=indata)
  rhsdata <- eval(rhs, envir=indata)

  # reformat
  if(is.vector(lhsdata)) {
    lhsdata <- matrix(lhsdata, ncol=1)
    colnames(lhsdata) <- paste(deparse(lhs), collapse="")
  }
  # check lhs names exist
  lnames <- colnames(lhsdata)
  nc <- ncol(lhsdata)
  lnames0 <- paste("V", seq_len(nc), sep="")
  if(length(lnames) != nc)
    colnames(lhsdata) <- lnames0
  else if(any(uhoh <- !nzchar(lnames)))
    colnames(lhsdata)[uhoh] <- lnames0[uhoh]

  # check rhs data
  if(is.matrix(rhsdata))
    stop("rhs of formula should yield a vector")
  rhsdata <- as.numeric(rhsdata)

  nplots <- ncol(lhsdata)
  allind <- 1:nplots

  # extra plots may be implied by 'shade'
  explicit.lhs.names <- colnames(lhsdata)
  
  if(!is.null(shade)) {
    # select columns by name or number
    names(allind) <- explicit.lhs.names
    shind <- try(allind[shade])
    if(inherits(shind, "try-error")) 
      stop("The argument shade should be a valid subset index for columns of x")
    if(any(nbg <- is.na(shind))) {
      # columns not included in formula; get them
      morelhs <- try(as.matrix(indata[ , shade[nbg], drop=FALSE]))
      if(inherits(morelhs, "try-error")) 
        stop("The argument shade should be a valid subset index for columns of x")
      nmore <- ncol(morelhs)
      lhsdata <- cbind(lhsdata, morelhs)
      shind[nbg] <- nplots + seq_len(nmore)
      nplots <- nplots + nmore
      lty <- c(lty, rep(lty[1], nmore))
      col <- c(col, rep(col[1], nmore))
      lwd <- c(lwd, rep(lwd[1], nmore))
    }
  }
  
  # restrict data to subset if desired
  if(!is.null(subset)) {
    keep <- if(is.character(subset))
		eval(parse(text=subset), envir=indata)
            else
                eval(subset, envir=indata)
    lhsdata <- lhsdata[keep, , drop=FALSE]
    rhsdata <- rhsdata[keep]
  } 

  # -------------------- determine plotting limits ----------------------
  
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

  # -------------  work out how to label the plot --------------------

  # extract plot labels 
  labl <- attr(x, "labl")
  # expand plot labels
  if(!is.null(fname <- attr(x, "fname")))
    labl <- sprintf(labl, fname)
  # create plot label map (key -> algebraic expression)
  map <- fvlabelmap(x) 

  # ......... label for x axis ..................

  if(is.null(xlab)) {
    argname <- fvnames(x, ".x")
    if(as.character(fmla)[3] == argname) {
      # the x axis is the default function argument.
      # Add name of unit of length 
      ax <- summary(unitname(x))$axis
      xlab <- if(!is.null(ax)) paste(argname, ax) else as.expression(as.name(argname)) 
    } else {
      # map ident to label
      xlab <- eval(substitute(substitute(rh, mp), list(rh=rhs, mp=map)))
    }
  }
  if(is.language(xlab) && !is.expression(xlab))
    xlab <- as.expression(xlab)
  
  # ......... label for y axis ...................

  leftside <- lhs.original
  if(ncol(lhsdata) > 1) {
    # For labelling purposes only, simplify the LHS by 
    # replacing 'cbind(.....)' by '.'
    # even if not all columns are included.
    leftside <- paste(as.expression(leftside))
    cb <- paste("cbind(",
                paste(explicit.lhs.names, collapse=", "),
                ")", sep="")
    compactleftside <- gsub(cb, ".", leftside, fixed=TRUE)
    # Separately expand "." to cbind(.....) and ".x", ".y" to their real names
    cball <- paste("cbind(",
                paste(fvnames(x, "."), collapse=", "),
                ")", sep="")
    expandleftside <- gsub(".x", fvnames(x, ".x"), leftside, fixed=TRUE)
    expandleftside <- gsub(".y", fvnames(x, ".y"), expandleftside, fixed=TRUE)
    expandleftside <- gsub(".", cball, expandleftside, fixed=TRUE)
    # convert back to language
    compactleftside <- as.formula(paste(compactleftside, "~1"))[[2]]
    expandleftside <- as.formula(paste(expandleftside, "~1"))[[2]]
  } else {
    compactleftside <- expandleftside <- leftside
  }

  # construct label for y axis
  if(is.null(ylab)) {
    yl <- attr(x, "yexp")
    if(defaultplot && !is.null(yl)) {
      ylab <- yl
    } else {
      # replace "." and short identifiers by plot labels
      ylab <- eval(substitute(substitute(le, mp),
                                list(le=compactleftside, mp=map)))
    }
  }
  if(is.language(ylab) && !is.expression(ylab))
    ylab <- as.expression(ylab)

  # ------------------ start plotting ---------------------------

  # create new plot
  if(!add)
    do.call("plot.default",
            resolve.defaults(list(xlim, ylim, type="n"),
                             list(xlab=xlab, ylab=ylab),
                             list(...),
                             list(main=xname)))


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

  if(!is.null(shade)) {
    # shade region between critical boundaries
    # extract relevant columns for shaded bands
    shdata <- lhsdata[, shind]
    if(!is.matrix(shdata) || ncol(shdata) != 2) 
      stop("The argument shade should select two columns of x")
    # determine plot limits for shaded bands
    shdata1 <- shdata[,1]
    shdata2 <- shdata[,2]
    rhsOK <- is.finite(rhsdata)
    shade1OK <- rhsOK & is.finite(shdata1)
    shade2OK <- rhsOK & is.finite(shdata2)
    shadeOK <- shade1OK & shade2OK
    # work out which one is the upper limit
    up1 <- all(shdata1[shadeOK] > shdata2[shadeOK])
    # half-infinite intervals
    if(!is.null(ylim)) {
      shdata1[shade2OK & !shade1OK] <- if(up1) ylim[2] else ylim[1]
      shdata2[shade1OK & !shade2OK] <- if(up1) ylim[1] else ylim[2]
      shadeOK <- shade1OK | shade2OK
    } 
    # plot grey polygon 
    polygon(c(rhsdata[shadeOK], rev(rhsdata[shadeOK])),
            c(shdata1[shadeOK],  rev(shdata2[shadeOK])),
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
  
  # ----------------- plot lines ------------------------------

  for(i in allind)
    lines(rhsdata, lhsdata[,i], lty=lty[i], col=col[i], lwd=lwd[i])

  # determine legend 
  if(nplots == 1)
    return(invisible(NULL))
  else {
    key <- colnames(lhsdata)
    mat <- match(key, names(x))
    keyok <- !is.na(mat)
    matok <- mat[keyok]
    legdesc <- rep("constructed variable", length(key))
    legdesc[keyok] <- attr(x, "desc")[matok]
    leglabl <- lnames0
    leglabl[keyok] <- labl[matok]
    ylab <- attr(x, "ylab")
    if(!is.null(ylab)) {
      if(is.language(ylab))
        ylab <- deparse(ylab)
      legdesc <- sprintf(legdesc, ylab)
    }
    # compute legend info
    legtxt <- key
    if(legendmath) {
      legtxt <- leglabl
      if(defaultplot) {
        # try to convert individual labels to expressions
        fancy <- try(parse(text=leglabl), silent=TRUE)
      } else {
        # try to navigate the parse tree
        fancy <- try(fvlegend(x, expandleftside), silent=TRUE)
      }
      if(!inherits(fancy, "try-error"))
        legtxt <- fancy
    }
    # plot legend
    if(!is.null(legend) && legend) 
      do.call("legend",
              resolve.defaults(legendargs,
                               list(x=legendpos, legend=legtxt, lty=lty, col=col),
                               list(inset=0.05, y.intersp=if(legendmath) 1.3 else 1),
                               .StripNull=TRUE))

    df <- data.frame(lty=lty, col=col, key=key, label=paste.expr(legtxt),
                      meaning=legdesc, row.names=key)
    return(df)
  }
}


