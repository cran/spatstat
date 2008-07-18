#
#	ppp.S
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.57 $	$Date: 2008/07/16 17:04:07 $
#
#	A point pattern contains the following entries:	
#
#		$window:	an object of class 'owin'
#				defining the observation window
#
#		$n:	the number of points (for efficiency)
#	
#		$x:	
#		$y:	vectors of length n giving the Cartesian
#			coordinates of the points.
#
#	It may also contain the entry:	
#
#		$marks:	a vector of length n
#			whose entries are interpreted as the
#			'marks' attached to the corresponding points.	
#	
#--------------------------------------------------------------------------
ppp <- function(x, y, ..., window, marks, check=TRUE ) {
  # Constructs an object of class 'ppp'
  #
  if(!missing(window))
    verifyclass(window, "owin")
  else
    window <- owin(...)
          
  n <- length(x)
  if(length(y) != n)
    stop("coordinate vectors x and y are not of equal length")
  
  # validate x, y coordinates
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  names(x) <- NULL
  names(y) <- NULL
  
  # check (x,y) points lie inside window
  if(check && n > 0) {
    ok <- inside.owin(x, y, window)
    nout <- sum(!ok)
    if(nout > 0) {
      warning(paste(nout,
                    ngettext(nout, "point was", "points were"),
                    "rejected as lying outside the specified window"))
      rr <- ripras(x,y)
      bb <- bounding.box.xy(x,y)
      bb <- bounding.box(rr, bb, window)
      rejectwindow <-
        if(!is.null(rr)) rebound.owin(rr, bb) else bb
      rejects <- ppp(x[!ok], y[!ok], window=rejectwindow, check=FALSE)
      # discard illegal points
      x <- x[ok]
      y <- y[ok]
      n <- length(x)
    }
  } else nout <- 0
  # initialise ppp object
  pp <- list(window=window, n=n, x=x, y=y)
  # add marks if any
  if(!missing(marks) && !is.null(marks)) {
    if(is.matrix(marks)) 
      marks <- as.data.frame(marks)
    if(is.data.frame(marks)) {
      pp$markformat <- "dataframe"
      if(.Spatstat.Forbids.Df)
        stop("Sorry, not implemented for data frame of marks",
             call.=FALSE)
      if(nout > 0) {
        marks <- marks[ok, ]
        marks(rejects) <- marks[!ok,]
      }
      if(nrow(marks) != n)
        stop("number of rows of marks != length of x and y")
    } else {
      pp$markformat <- "vector"
      if(nout > 0) {
        marks(rejects) <- marks[!ok]
        marks <- marks[ok]
      }
      if(length(marks) != n)
        stop("length of marks vector != length of x and y")
      names(marks) <- NULL
    }
    pp$marks <- marks
  } else {
    pp$markformat <- "none"
  }
  class(pp) <- "ppp"
  if(check && any(duplicated(pp)))
    warning("data contain duplicated points")
  if(nout > 0) 
    attr(pp, "rejects") <- rejects
  pp
}

#
#--------------------------------------------------------------------------
#

is.ppp <- function(x) { inherits(x, "ppp") }

#
#--------------------------------------------------------------------------
#

as.ppp <- function(X, ..., fatal=TRUE) {
  UseMethod("as.ppp")
}

as.ppp.ppp <- function(X, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=FALSE))$check
  return(ppp(X$x, X$y, window=X$window, marks=X$marks, check=check))
}

as.ppp.quad <- function(X, ..., fatal=TRUE) {
  return(union.quad(X))
}

as.ppp.data.frame <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(!(ncol(X) %in% 2:3)) {
    if(fatal)
      stop("X must be a two-column or three-column data frame")
    else
      return(NULL)
  }

  if(is.null(W)) {
    if(fatal)
      stop("x,y coords given but no window specified")
    else
      return(NULL)
  } 
    
  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, check=check)
  }

  # add marks from third column
  if(ncol(X) == 3)
    Z <- Z %mark% X[,3]

  return(Z)
}
    
as.ppp.matrix <- function(X, W = NULL, ..., fatal=TRUE) {
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(!verifyclass(X, "matrix", fatal=fatal)
     || !is.numeric(X)) {
    if(fatal)
      stop("X must be a two-column matrix of numbers")
    else
      return(NULL)
  }
  if(!(ncol(X) %in% 2:3)) {
    if(fatal)
      stop("X must be a two-column or three-column matrix of numbers")
    else
      return(NULL)
  }

  if(is.null(W)) {
    if(fatal)
      stop("x,y coords given but no window specified")
    else
      return(NULL)
  } 
    
  if(is.function(W))
    Z <- cobble.xy(X[,1], X[,2], W, fatal)
  else {
    win <- as.owin(W)
    Z <- ppp(X[,1], X[,2], window = win, check=check)
  }

  # add marks from third column
  if(ncol(X) == 3)
    Z <- Z %mark% X[,3]

  return(Z)
}
    
as.ppp.default <- function(X, W=NULL, ..., fatal=TRUE) {
	# tries to coerce data X to a point pattern
	# X may be:
	#	1. a structure with entries x, y, xl, xu, yl, yu
	#	2. a structure with entries x, y, area where
        #                    'area' has entries xl, xu, yl, yu
	#	3. a structure with entries x, y
	# The second argument W is coerced to an object of class 'owin' by the 
	# function "as.owin" in window.S
        # If X also has an entry X$marks
        # then this will be interpreted as the marks vector for the pattern.
	#
  check <- resolve.defaults(list(...), list(check=TRUE))$check
  if(checkfields(X, c("x", "y", "xl", "xu", "yl", "yu"))) {
		xrange <- c(X$xl, X$xu)
		yrange <- c(X$yl, X$yu)
		if(is.null(X$marks))
			Z <- ppp(X$x, X$y, xrange, yrange, check=check)
		else
			Z <- ppp(X$x, X$y, xrange, yrange, 
				marks=X$marks, check=check)
		return(Z)
        } else if(checkfields(X, c("x", "y", "area"))
                  && checkfields(X$area, c("xl", "xu", "yl", "yu"))) {
                win <- as.owin(X$area)
                if (is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks = X$marks, check=check)
                return(Z)
	} else if(checkfields(X, c("x", "y"))) {
                if(is.function(W))
                  return(cobble.xy(X$x, X$y, W, fatal))
		if(is.null(W)) {
                  if(fatal)
                    stop("x,y coords given but no window specified")
                  else
                    return(NULL)
                }
		win <- as.owin(W)
		if(is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win, check=check)
                else
                  Z <- ppp(X$x, X$y, window=win, marks=X$marks, check=check)
                return(Z)
	} else {
          if(fatal)
            stop("Can't interpret X as a point pattern")
          else
            return(NULL)
        }
}

cobble.xy <- function(x, y, f=ripras, fatal=TRUE) {
  if(!is.function(f))
    stop("f is not a function")
  w <- f(x,y)
  if(!is.owin(w)) {
    gripe <- "Supplied function f did not return an owin object"
    if(fatal)
      stop(gripe)
    else {
      warning(gripe)
      return(NULL)
    }
  }
  return(ppp(x, y, window=w))
}
  

# --------------------------------------------------------------

"[.ppp" <-
"subset.ppp" <-
  function(x, i, j, drop, ...) {

        verifyclass(x, "ppp")
        
        if(missing(i) && missing(j))
          return(x)

        if(!missing(i)) {
          if(inherits(i, "owin")) {
            # i is a window
            window <- i
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else if(inherits(i, "im")) {
            # i is an image
            if(i$type != "logical")
              stop(paste("Subset operator X[i] undefined",
                         "when i is a pixel image",
                         "unless it has logical values"), call.=FALSE)
            # convert logical image to window
            e <- sys.frame(sys.nframe())
            window <- solutionset(i, e)
            ok <- inside.owin(x$x, x$y, window)
            x <- ppp(x$x[ok], x$y[ok], window=window, #SIC
                     marks=marksubset(x$marks, ok),
                     check=FALSE)
          } else {
            # assume i is a subset index
            if(x$n == 0)
              return(x)
            subset <- i
            x <- ppp(x$x[subset], x$y[subset], window=x$window,
                     marks=marksubset(x$marks, subset),
                     check=FALSE)
          } 
        }

        if(!missing(j))
          x <- x[j]   # invokes code above

        return(x)
}


# ------------------------------------------------------------------

cut.ppp <- function(x, ...) {
  x <- as.ppp(x)
  if(!is.marked(x))
    stop("x has no marks to cut")
  marks(x) <- cut(marks(x, dfok=FALSE), ...)
  return(x)
}

# ------------------------------------------------------------------
#
#
scanpp <- function(filename, window, header=TRUE, dir="", multitype=FALSE) {
  filename <- if(dir=="") filename else
              paste(dir, filename, sep=.Platform$file.sep)
  df <- read.table(filename, header=header)
  if(header) {
    x <- df$x
    y <- df$y
    colnames <- dimnames(df)[[2]]
    xycolumns <- match(colnames, c("x","y"), 0)
  } else {
    # assume x, y given in columns 1, 2 respectively
    x <- df[,1]
    y <- df[,2]
    xycolumns <- c(1,2)
  }
  if(ncol(df) == 2) 
      X <- ppp(x, y, window=window)
  else {
    marks <- df[ , -xycolumns]
    if(multitype) 
      marks <- factor(marks)
    X <- ppp(x, y, window=window, marks = marks)
  }
  X
}

#-------------------------------------------------------------------

"markspace.integral" <-
  function(X) {
  verifyclass(X, "ppp")
  if(!is.marked(X))
    return(1)
  if(is.multitype(X))
    return(length(levels(marks(X))))
  else
    stop("Don't know how to compute total mass of mark space")
}

#-------------------------------------------------------------------

print.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  ism <- is.marked(x, dfok=TRUE)
  cat(paste(if(ism) "marked" else NULL,
            "planar point pattern:",
            x$n,
            ngettext(x$n, "point", "points"),
            "\n"))
  if(ism) {
    mks <- marks(x, dfok=TRUE)
    if(is.data.frame(mks)) {
      # data frame of marks
      cat(paste("Mark variables: ", paste(names(mks), collapse=", "), "\n"))
    } else {
      # vector of marks
      if(is.factor(mks)) {
        cat("multitype, with ")
        cat(paste("levels =", paste(levels(mks), collapse="\t"),"\n"))
      } else {
        cat(paste("marks are",
                  if(is.numeric(mks)) "numeric,",
                  "of type", sQuote(typeof(mks)), "\n"))
      }
    }
  }
  print(x$window)
  if(!is.null(rejects <- attr(x, "rejects"))) {
    nrejects <- rejects$n
    cat(paste("\n***",
              nrejects,
              ngettext(nrejects, "illegal point", "illegal points"),
              "stored in",
              paste("attr(,", dQuote("rejects"), ")", sep=""),
              "***\n"))
  }
  return(invisible(NULL))
}

summary.ppp <- function(object, ..., checkdup=TRUE) {
  verifyclass(object, "ppp")
  result <- list()
  result$is.marked <- is.marked(object, dfok=TRUE)
  result$n <- object$n
  result$window <- summary(object$window)
  result$intensity <- result$n/result$window$area
  if(checkdup)
    result$nduplicated <- sum(duplicated(object))
  if(result$is.marked) {
    mks <- marks(object, dfok=TRUE)
    if(result$multiple.marks <- is.data.frame(mks)) {
      result$marknames <- names(mks)
      result$is.numeric <- FALSE
      result$marktype <- "dataframe"
      result$is.multitype <- FALSE
    } else {
      result$is.numeric <- is.numeric(mks)
      result$marknames <- "marks"
      result$marktype <- typeof(mks)
      result$is.multitype <- is.multitype(object)
    }
    if(result$is.multitype) {
      tm <- as.vector(table(mks))
      tfp <- data.frame(frequency=tm,
                        proportion=tm/sum(tm),
                        intensity=tm/result$window$area,
                        row.names=levels(mks))
      result$marks <- tfp
    } else 
      result$marks <- summary(mks)
  }
  class(result) <- "summary.ppp"
  if(!is.null(rejects <- attr(object, "rejects"))) 
    result$rejects <- rejects$n
  return(result)
}

print.summary.ppp <- function(x, ..., dp=3) {
  verifyclass(x, "summary.ppp")
  cat(paste(if(x$is.marked) "Marked planar " else "Planar ",
            "point pattern: ",
            x$n,
            " points\n",
            sep=""))
  oneline <- resolve.defaults(list(...), list(oneline=FALSE))$oneline
  if(oneline) return(invisible(NULL))
  unitinfo <- summary(x$window$units)
  cat(paste("Average intensity",
            signif(x$intensity,dp),
            "points per square",
            unitinfo$singular,
            unitinfo$explain,
            "\n"))
  ndup <- x$nduplicated
  if((!is.null(ndup)) && (ndup > 0))
    cat("\n*Pattern contains duplicated points*\n")
  if(x$is.marked) {
    if(x$multiple.marks) {
      cat(paste("Mark variables: ", paste(x$marknames, collapse=", "), "\n"))
      cat("Summary:\n")
      print(x$marks)
    } else if(x$is.multitype) {
      cat("Multitype:\n")
      print(signif(x$marks,dp))
    } else {
      cat(paste("marks are ",
                if(x$is.numeric) "numeric, ",
                "of type", sQuote(x$marktype), "\n"))
      cat("Summary:\n")
      print(x$marks)
    }
  }
  cat("\n")
  print(x$window)
  if(!is.null(nrejects <- x$rejects)) 
    cat(paste("\n***",
              nrejects,
              ngettext(nrejects, "illegal point", "illegal points"),
              "stored in",
              paste("attr(,", dQuote("rejects"), ")", sep=""),
              "***\n"))
  return(invisible(x))
}

# ---------------------------------------------------------------

"%mark%" <-
setmarks <- function(X, m) {
  verifyclass(X, "ppp")
  if(is.null(m))
    return(unmark(X))
  if(!is.data.frame(m)) {
    # vector of marks
    if(length(m) == 1) m <- rep(m, X$n)
    else if(X$n == 0) m <- rep(m, 0) # ensures marked pattern is obtained
    else if(length(m) != X$n) stop("number of points != number of marks")
    Y <- ppp(X$x,X$y,window=X$window,marks=m, check=FALSE)
    return(Y)
  }
  # data frame of marks
  if(.Spatstat.Forbids.Df) 
    stop("Sorry, not implemented for data frame of marks", call.=FALSE)
  if(nrow(m) == 1 || X$n == 0) {
    # replicate data frame
    m <- as.data.frame(lapply(as.list(m),
                              function(x, k) { rep(x, k) },
                              k=X$n))
  } else if(nrow(m) != X$n)
    stop("number of rows of data frame != number of points")
  Y <- ppp(X$x,X$y,window=X$window,marks=m, check=FALSE)
}

identify.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  if(!is.marked(x) || "labels" %in% names(list(...)))
    identify(x$x, x$y, ...)
  else {
    marx <- marks(x, dfok=FALSE)
    marques <- if(is.numeric(marx)) paste(signif(marx, 3)) else paste(marx)
    id <- identify(x$x, x$y, labels=marques, ...)
    mk <- marx[id]
    if(is.factor(marx)) mk <- levels(marx)[mk]
    cbind(id=id, marks=mk)
  }
}

rebound.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  x$window <- rebound.owin(x$window, ...)
  return(x)
}

as.data.frame.ppp <- function(x, row.names=NULL, ...) {
  if(!is.marked(x))
    data.frame(x=x$x, y=x$y, row.names=row.names)
  else
    data.frame(x=x$x, y=x$y, marks=x$marks, row.names=row.names)
}

