#
#	ppp.S
#
#	A class 'ppp' to define point patterns
#	observed in arbitrary windows in two dimensions.
#
#	$Revision: 4.6 $	$Date: 2002/05/13 12:41:10 $
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

as.ppp <- function(X, W = NULL) {
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
		if(is.null(W))
			stop("x,y coords given but no window specified")
		win <- as.owin(W)
		Z <- ppp(X[,1], X[,2], window = win)
		return(Z)
	} else if(checkfields(X, c("x", "y"))) {
		if(is.null(W))
			stop("x,y coords given but no window specified")
		win <- as.owin(W)
		if(is.null(X$marks))
                  Z <- ppp(X$x, X$y, window=win)
                else
                  Z <- ppp(X$x, X$y, window=win, marks=X$marks)
                return(Z)
	} else
		stop("Can't interpret X as a point pattern")
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


"is.marked.ppp" <-
function(X, na.action="warn") {
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
