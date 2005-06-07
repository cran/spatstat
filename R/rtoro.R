#
#   rtoro.R
#
#   random toroidal shift
#
#   $Revision: 1.1 $   $Date: 2005/06/08 01:33:39 $
#
#
rtoro <- function(X, which=NULL) {

  if(inherits(X, "splitppp")) {
    if(is.null(which))
      which <- seq(X)
    Xsub <- X[which]
    shiftXsub <- lapply(Xsub, rtoro)
    X[which] <- shiftXsub
    return(X)
  }

  verifyclass(X, "ppp")

  if(!is.null(which)) {
    if(!is.marked(X) || !is.factor(X$marks))
      stop("pattern is not multitype: use of \`which\' undefined")
    split(X) <- rtoro(split(X), which)
    return(X)
  }

  # vanilla point pattern
  # shift all points
  
  W <- X$window
  W <- rescue.rectangle(W)
  if(W$type != "rectangle")
    stop("Only meaningful for rectangular windows")

  xr <- W$xrange
  yr <- W$yrange
  width <- diff(xr)
  height <- diff(yr)
  
  xjump <- runif(1, min=0, max=width)
  yjump <- runif(1, min=0, max=height)

  x <- X$x + xjump
  y <- X$y + yjump

  x <- ifelse(x < xr[1], x + width,  ifelse(x > xr[2], x - width, x))
  y <- ifelse(y < yr[1], y + height, ifelse(y > yr[2], y - height, y))

  X$x <- x
  X$y <- y
  
  return(X)
}

