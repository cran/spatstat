#
#	ppp.S
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.23 $	$Date: 2005/03/09 19:49:55 $
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
ppp <- function(x, y, ..., window, marks ) {
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
        # initialise ppp object
	pp <- list(window=window, n=n, x=x, y=y)
        # add marks if any
	if(!missing(marks) && !is.null(marks)) {
                if(is.matrix(marks) || is.data.frame(marks))
                  stop(paste("Attempted to create point pattern with",
                             ncol(marks), "columns of mark data;",
                             "multidimensional marks are not yet implemented"))
		if(length(marks) != n)
			stop("length of marks vector != length of x and y")
                names(marks) <- NULL
		pp$marks <- marks
	}
	class(pp) <- "ppp"
	pp
}
#
#--------------------------------------------------------------------------
#

is.ppp <- function(x) { inherits(x, "ppp") }

as.ppp <- function(X, W = NULL, fatal=TRUE) {
	# tries to coerce data X to a point pattern
	# X may be:
	#	1. an object of class 'ppp'
	#	2a. a structure with entries x, y, xl, xu, yl, yu
	#	2b. a structure with entries x, y, area where
        #                    'area' has entries xl, xu, yl, yu
	#	3. a two-column matrix
	#	4. a structure with entries x, y
        #       5. a quadrature scheme (object of class 'quad')
	# In cases 3 and 4, we need the second argument W
	# which is coerced to an object of class 'owin' by the 
	# function "as.owin" in window.S
        # In cases 2 and 4, if X also has an entry X$marks
        # then this will be interpreted as the marks vector for the pattern.
	#
	if(verifyclass(X, "ppp", fatal=FALSE))
		return(X)
        else if(verifyclass(X, "quad", fatal=FALSE))
                return(union.quad(X))
	else if(checkfields(X, 	c("x", "y", "xl", "xu", "yl", "yu"))) {
		xrange <- c(X$xl, X$xu)
		yrange <- c(X$yl, X$yu)
		if(is.null(X$marks))
			Z <- ppp(X$x, X$y, xrange, yrange)
		else
			Z <- ppp(X$x, X$y, xrange, yrange, 
				marks=X$marks)
		return(Z)
        } else if(checkfields(X, c("x", "y", "area"))
                  && checkfields(X$area, c("xl", "xu", "yl", "yu"))) {
                win <- as.owin(X$area)
                if (is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win)
                else
                  Z <- ppp(X$x, X$y, window=win, marks = X$marks)
                return(Z)
	} else if(is.matrix(X) && is.numeric(X) && ncol(X) == 2) {
                if(is.function(W))
                  return(cobble.ppp(X[,1], X[,2], W, fatal))
		if(is.null(W)) {
                  if(fatal)
                    stop("x,y coords given but no window specified")
                  else
                    return(NULL)
                } 
                win <- as.owin(W)
		Z <- ppp(X[,1], X[,2], window = win)
		return(Z)
	} else if(checkfields(X, c("x", "y"))) {
                if(is.function(W))
                  return(cobble.ppp(X$x, X$y, W, fatal))
		if(is.null(W)) {
                  if(fatal)
                    stop("x,y coords given but no window specified")
                  else
                    return(NULL)
                }
		win <- as.owin(W)
		if(is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win)
                else
                  Z <- ppp(X$x, X$y, window=win, marks=X$marks)
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
  function(x, subset, window, drop, ...) {

        verifyclass(x, "ppp")

        trim <- !missing(window)
        thin <- !missing(subset)
        if(!thin && !trim)
          stop("Please specify a subset (to thin the pattern) or a window (to trim it)")

        # thin first, according to 'subset'
        if(!thin)
          Y <- x
        else
          Y <- ppp(x$x[subset],
                   x$y[subset],
                   window=x$window,
                   marks=if(is.null(x$marks)) NULL else x$marks[subset])

        # now trim to window 
        if(trim) {
          ok <- inside.owin(Y$x, Y$y, window)
          Y <- ppp(Y$x[ok], Y$y[ok],
                   window=window,  # SIC
                   marks=if(is.null(Y$marks)) NULL else Y$marks[ok])
        }
        
        return(Y)
}

# ------------------------------------------------------------------

cut.ppp <- function(x, ...) {
  x <- as.ppp(x)
  if(!is.marked(x))
    stop("x has no marks to cut")
  x$marks <- cut(x$marks, ...)
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

"is.marked.ppp" <-
function(X, na.action="warn", ...) {
    verifyclass(X, "ppp")
    if(is.null(X$marks))
      return(FALSE)
    if(any(is.na(X$marks)))
      switch(na.action,
             warn = {
               warning(paste("some mark values are NA in the point pattern",
                    deparse(substitute(X))))
             },
             fatal = {
               return(FALSE)
             },
             ignore = {
               return(TRUE)
             }
      )
    return(TRUE)
}

"is.marked" <-
function(X, ...) {
  UseMethod("is.marked")
}

"is.marked.default" <-
  function(...) { return(FALSE) }

"unmark" <-
function(X) {
  if(inherits(X, "ppp")) {
    X$marks <- NULL
    return(X)
  } else if(inherits(X, "splitppp")) {
    Y <- lapply(X, unmark)
    class(Y) <- c("splitppp", class(Y))
    return(Y)
  } else
    stop(paste("X must be a point pattern (class \"ppp\")\n",
               "or a list of point patterns (class \"splitppp\")"))
}

"markspace.integral" <-
  function(X) {
  verifyclass(X, "ppp")
  if(!is.marked(X))
    return(1)
  if(is.factor(X$marks))
    return(length(levels(X$marks)))
  else
    stop("Don't know how to compute total mass of mark space")
}

#-------------------------------------------------------------------

print.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  ism <- is.marked(x)
  cat(paste(if(ism) "marked" else NULL,
            "planar point pattern:",
            x$n,
            "points\n"))
  if(ism) {
    mks <- x$marks
    if(is.factor(mks)) {
      cat("multitype, with ")
      cat(paste("levels =", paste(levels(mks), collapse="\t"),"\n"))
    } else {
      cat(paste("marks are ",
                if(is.numeric(mks)) "numeric, ",
                "of type \"", typeof(mks),
                "\"\n", sep=""))
    }
  }
  print(x$window)
  return(invisible(NULL))
}

summary.ppp <- function(object, ...) {
  verifyclass(object, "ppp")
  result <- list()
  result$is.marked <- is.marked(object)
  result$n <- object$n
  result$window <- summary(object$window)
  result$intensity <- result$n/result$window$area
  if(result$is.marked) {
    mks <- object$marks
    result$is.multitype <- is.factor(mks)
    result$is.numeric <- is.numeric(mks)
    result$marktype <- typeof(mks)
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
  return(result)
}

print.summary.ppp <- function(x, ..., dp=3) {
  verifyclass(x, "summary.ppp")
  cat(paste(if(x$is.marked) "Marked planar " else "Planar ",
            "point pattern: ",
            x$n,
            " points\n",
            sep=""))
  cat(paste("Average intensity",
            signif(x$intensity,dp), "points per unit area\n"))
  if(x$is.marked)
    if(x$is.multitype) {
      cat("Marks:\n")
      print(signif(x$marks,dp))
    } else {
      cat(paste("marks are ",
                if(x$is.numeric) "numeric, ",
                "of type \"", x$marktype,
                "\"\n", sep=""))
      cat("summary:\n")
      print(x$marks)
    }
  cat("\n")
  print(x$window)
  return(invisible(x))
}

# ---------------------------------------------------------------

"%mark%" <-
setmarks <- function(X, m) {
  verifyclass(X, "ppp")
  if(length(m) == 1) m <- rep(m, X$n)
  else if(X$n == 0) m <- rep(m, 0) # ensures marked pattern is obtained
  else if(length(m) != X$n) stop("number of points != number of marks")
  Y <- ppp(X$x,X$y,window=X$window,marks=m)
  return(Y)
}

identify.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  if(!is.marked(x) || "labels" %in% names(list(...)))
    identify(x$x, x$y, ...)
  else {
    marx <- x$marks
    marques <- if(is.numeric(marx)) paste(signif(marx, 3)) else paste(marx)
    id <- identify(x$x, x$y, labels=marques, ...)
    mk <- marx[id]
    if(is.factor(marx)) mk <- levels(marx)[mk]
    cbind(id=id, marks=mk)
  }
}


