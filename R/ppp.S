#
#	ppp.S
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.12 $	$Date: 2003/05/02 10:55:09 $
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
	pp <- list(window=window, n=n, x=x, y=y)
	if(!missing(marks) && !is.null(marks)) {
                if(is.matrix(marks) || is.data.frame(marks))
                  stop(paste("Attempted to create point pattern with",
                             ncol(marks), "columns of mark data;",
                             "multidimensional marks are not yet implemented"))
		if(length(marks) != n)
			stop("length of marks vector != length of x and y")
		pp$marks <- marks
	}
	class(pp) <- "ppp"
	pp
}
#
#--------------------------------------------------------------------------
#

as.ppp <- function(X, W = NULL, fatal=TRUE) {
	# tries to coerce data X to a point pattern
	# X may be:
	#	1. an object of class 'ppp'
	#	2. a structure with entries x, y, xl, xu, yl, yu
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
	} else if(is.matrix(X) && is.numeric(X)) {
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
  filename <- paste(dir, filename, sep="")
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


"superimpose" <-
  function(...)
{
  # superimpose any number of point patterns
  # ASSUMED TO BE IN THE SAME WINDOW
  # WITH THE SAME MARK SPACE if relevant
  
  # concatenate lists of (x,y) coordinates
  XY <- concatxy(...)
  # concatenate vectors of marks
  M <- unlist(lapply(list(...), function(x) {x$marks}))
  if(length(M) > 0 && length(M) != length(XY$x)) {
    warning("marks not present in all patterns -- ignored them.")
    M <- NULL
  }
  # determine window
  P <- ..1
  verifyclass(P, "ppp")
  win <- P$window
  #
  # determine type of marks
  if(!is.null(M) && is.factor(P$marks)) {
    M <- factor(M)
    levels(M) <- levels(P$marks)
  }
  #
  ppp(XY$x, XY$y, window=win, marks=M)
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
  verifyclass(X, "ppp")
  X$marks <- NULL
  X
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
  ism <- is.marked(object)
  cat(paste(if(ism) "marked " else NULL,
            "planar point pattern: ",
            object$n,
            " points\n",
            sep=""))
  a <- area.owin(object$window)
  cat(paste("average intensity",
            signif(object$n/a,3), "points per unit area\n"))
  if(ism) {
    mks <- object$marks
    if(is.factor(mks)) {
      cat("Marks:\n")
      tm <- as.vector(table(mks))
      tfp <- data.frame(frequency=tm,
                        proportion=signif(tm/sum(tm),3),
                        intensity=signif(tm/a, 3),
                        row.names=levels(mks))
      print(tfp)
    } else {
      cat(paste("marks are ",
                if(is.numeric(mks)) "numeric, ",
                "of type \"", typeof(mks),
                "\"\n", sep=""))
      cat("summary:\n")
      print(summary(mks))
    }
  }
  summary(object$window)
  return(invisible(NULL))
}

# ---------------------------------------------------------------

"%mark%" <-
setmarks <- function(X, m) {
  verifyclass(X, "ppp")
  Y <- ppp(X$x,X$y,window=X$window,marks=m)
  return(Y)
}

identify.ppp <- function(x, ...) {
  verifyclass(x, "ppp")
  if(!is.marked(x) || "labels" %in% names(list(...)))
    identify(x$x, x$y, ...)
  else {
    id <- identify(x$x, x$y, labels=paste(x$marks), ...)
    mk <- x$marks[id]
    if(is.factor(x$marks)) mk <- levels(x$marks)[mk]
    cbind(id=id, marks=mk)
  }
}


