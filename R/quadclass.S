#
#	quadclass.S
#
#	Class 'quad' to define quadrature schemes
#	in (rectangular) windows in two dimensions.
#
#	$Revision: 4.2 $	$Date: 2002/01/17 09:44:54 $
#
# An object of class 'quad' contains the following entries:
#
#	$data:	an object of class 'ppp'
#		defining the OBSERVATION window, 
#		giving the locations (& marks) of the data points.
#
#	$dummy:	object of class 'ppp'
#		defining the QUADRATURE window, 
#		giving the locations (& marks) of the dummy points.
#	
#	$w: 	vector giving the nonnegative weights for the
#		data and dummy points (data first, followed by dummy)
#
#		w may also have an attribute attr(w, "zeroes")
#               equivalent to (w == 0). If this is absent
#               then all points are known to have positive weights.
#
#       The combined (data+dummy) vectors of x, y coordinates of the points, 
#       and their weights, are extracted using standard functions 
#       x.quad(), y.quad(), w.quad() etc.
#
#-------------------------------------------------------------

quad <- function(data, dummy, w) {
  
  data <- as.ppp(data)
  dummy <- as.ppp(dummy)

  n <- data$n + dummy$n
	
  if(missing(w))
    w <- rep(1, n)
  else if(length(w) != n)
    stop("length of weights vector w is not equal to total number of points")

  if(is.null(attr(w, "zeroes")) && any( w == 0))
	attr(w, "zeroes") <- (w == 0)

  Q <- list(data=data, dummy=dummy, w=w)
  class(Q) <- "quad"

  invisible(Q)
}

# ------------------ extractor functions ----------------------

x.quad <- function(Q) {
  verifyclass(Q, "quad")
  c(Q$data$x, Q$dummy$x)
}

y.quad <- function(Q) {
  verifyclass(Q, "quad")
  c(Q$data$y, Q$dummy$y)
}

w.quad <- function(Q) {
  verifyclass(Q, "quad")
  Q$w
}

n.quad <- function(Q) {
  verifyclass(Q, "quad")
  Q$data$n + Q$dummy$n
}

marks.quad <- function(Q) {
  verifyclass(Q, "quad")
  mdat <- Q$data$marks
  mdum <- Q$dummy$marks
  if(is.null(mdat) && is.null(mdum))
    return(NULL)
  if(is.null(mdat))
    mdat <- rep(NA, Q$data$n)
  if(is.null(mdum))
    mdum <- rep(NA, Q$dummy$n)
  mall <- c(mdat, mdum)
  if(is.factor(mdat) && is.factor(mdum) && all(levels(mdat) == levels(mdum))) {
    mall <- factor(mall)
    levels(mall) <- levels(mdat)
  }
  return(mall)
}

is.data <- function(Q) {
  verifyclass(Q, "quad")
  return(c(rep(T, Q$data$n),
	   rep(F, Q$dummy$n)))
}

equals.quad <- function(Q) {
    # return matrix E such that E[i,j] = (X[i] == U[j])
    # where X = Q$data and U = union.quad(Q)
    n <- Q$data$n
    m <- Q$dummy$n
    E <- matrix(F, nrow=n, ncol=n+m)
    diag(E) <- T
    E
}

  
union.quad <- function(Q) {
  verifyclass(Q, "quad")
  ppp(x= c(Q$data$x, Q$dummy$x),
      y= c(Q$data$y, Q$dummy$y),
      window=Q$dummy$window,
      marks=marks.quad(Q))
}
	
