#
#   rtoro.R
#
#   random toroidal shift
#
#   $Revision: 1.1 $   $Date: 2005/06/08 01:33:39 $
#
#
rtoro <- function(X, which=NULL, radius=NULL, width=NULL, height=NULL) {

  if(!is.null(radius) && !(is.null(width) && is.null(height)))
    stop("\`radius\' is incompatible with \`width\' and \`height\'")

  if(inherits(X, "splitppp")) {
    if(is.null(which))
      which <- seq(X)
    Xsub <- X[which]
    shiftXsub <- lapply(Xsub, rtoro, radius=radius, width=width, height=height)
    X[which] <- shiftXsub
    return(X)
  }

  verifyclass(X, "ppp")

  if(!is.null(which)) {
    if(!is.marked(X) || !is.factor(X$marks))
      stop("pattern is not multitype: use of \`which\' undefined")
    split(X) <- rtoro(split(X), which,
                      radius=radius, width=width, height=height)
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
  Wide <- diff(xr)
  High <- diff(yr)

  # generate random translation vector
  
  if(!is.null(radius)) 
    jump <- runifdisc(1, r=radius)
  else {
    if(is.null(width)) width <- Wide
    if(is.null(height)) height <- High
    jump <- list(x=runif(1, min=0, max=width),
                 y=runif(1, min=0, max=height))
  }

  # translate points
  x <- X$x + jump$x
  y <- X$y + jump$y

  # wrap points
  x <- ifelse(x < xr[1], x + Wide,  ifelse(x > xr[2], x - Wide, x))
  y <- ifelse(y < yr[1], y + High, ifelse(y > yr[2], y - High, y))

  # pack and ship
  X$x <- x
  X$y <- y
  
  return(X)
}

