#
#    xypolygon.S
#
#    $Revision: 1.3 $    $Date: 2002/01/17 10:01:49 $
#
#    low-level functions defined for polygons in list(x,y) format
#
verify.xypolygon <- function(p, fatal=T) {
  whinge <- NULL
  if(!is.list(p) || length(p) != 2 ||
     length(names(p)) != 2 || any(sort(names(p)) != c("x","y")))
    whinge <- "polygon must be a list with two components x and y"
  else if(is.null(p$x) || is.null(p$y) || !is.numeric(p$x) || !is.numeric(p$y))
    whinge <- "components x and y must be numeric vectors"
  else if(length(p$x) != length(p$y))
    whinge <- "lengths of x and y vectors unequal"
  ok <- is.null(whinge)
  if(!ok && fatal)
    stop(whinge)
  return(ok)
}

inside.xypolygon <- function(pts, polly, test01=T) {
  # pts:  list(x,y) points to be tested
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  # test01: logical - if TRUE, test whether all values in output are 0 or 1
  verify.xypolygon(pts)
  verify.xypolygon(polly)
  
  x <- pts$x
  y <- pts$y
  xp <- polly$x
  yp <- polly$y

  npts <- length(x)
  nedges <- length(xp)   # sic

  score <- rep(0, npts)
  on.boundary <- rep(F, npts)

  for(i in 1:nedges) {
    x0 <- xp[i]
    y0 <- yp[i]
    x1 <- if(i == nedges) xp[1] else xp[i+1]
    y1 <- if(i == nedges) yp[1] else yp[i+1]
    dx <- x1 - x0
    dy <- y1 - y0
    if(dx < 0) {
      # upper edge
      xcriterion <- (x - x0) * (x - x1)
      consider <- (xcriterion <= 0)
      if(any(consider)) {
        ycriterion <- y[consider] * dx - x[consider] * dy + (x0 * dy - y0 * dx)
        # closed inequality
        contrib <- (ycriterion >= 0) * ifelse(xcriterion[consider] == 0, 1/2, 1)
        # positive edge sign
        score[consider] <- score[consider] + contrib
        # detect whether any point lies on this segment
        on.boundary[consider] <- on.boundary[consider] | (ycriterion == 0)
      }
    } else if(dx > 0) {
      # lower edge
      xcriterion <- (x - x0) * (x - x1)
      consider <- (xcriterion <= 0)
      if(any(consider)) {
        ycriterion <- y[consider] * dx - x[consider] * dy + (x0 * dy - y0 * dx)
        # open inequality
        contrib <- (ycriterion < 0) * ifelse(xcriterion[consider] == 0, 1/2, 1)
        # negative edge sign
        score[consider] <- score[consider] - contrib
        # detect whether any point lies on this segment
        on.boundary[consider] <- on.boundary[consider] | (ycriterion == 0)
      }
    } else {
      # vertical edge
      consider <- (x == x0)
      if(any(consider)) {
        # zero score
        # detect whether any point lies on this segment
        yconsider <- y[consider]
        ycriterion <- (yconsider - y0) * (yconsider - y1)
        on.boundary[consider] <- on.boundary[consider] | (ycriterion <= 0)
      }
    }
  }
  
  # any point recognised as lying on the boundary gets score 1.
  score[on.boundary] <- 1

  if(test01) {
    # check sanity
    if(!all((score == 0) | (score == 1)))
      warning("internal error: some scores are not equal to 0 or 1")
  }

  attr(score, "on.boundary") <- on.boundary
  
  return(score)
}

area.xypolygon <- function(polly) {
  #
  # polly: list(x,y) vertices of a single polygon (n joins to 1)
  #
  verify.xypolygon(polly)
  
  xp <- polly$x
  yp <- polly$y
  
  nedges <- length(xp)   # sic
  
  # place x axis below polygon
  yp <- yp - min(yp) 

  # join vertex n to vertex 1
  nxt <- c(2:nedges, 1)

  # x step, WITH sign
  dx <- xp[nxt] - xp

  # average height 
  ym <- (yp + yp[nxt])/2
  
  -sum(dx * ym)
}

bdrylength.xypolygon <- function(polly) {
  verify.xypolygon(polly)
  xp <- polly$x
  yp <- polly$y
  nedges <- length(xp)
  nxt <- c(2:nedges, 1)
  dx <- xp[nxt] - xp
  dy <- yp[nxt] - yp
  sum(sqrt(dx^2 + dy^2))
}

reverse.xypolygon <- function(p) {
  # reverse the order of vertices
  # (=> change sign of polygon)

  verify.xypolygon(p)
  
  n <- length(p$x)

  return(list(x=p$x[n:1], y=p$y[n:1]))
}

overlap.xypolygon <- function(P, Q) {
  # compute area of overlap of two simple closed polygons 
  verify.xypolygon(P)
  verify.xypolygon(Q)
  
  xp <- P$x
  yp <- P$y
  np <- length(xp)
  nextp <- c(2:np, 1)

  xq <- Q$x
  yq <- Q$y
  nq <- length(xq)
  nextq <- c(2:nq, 1)

  # adjust y coordinates so all are nonnegative
  ylow <- min(c(yp,yq))
  yp <- yp - ylow
  yq <- yq - ylow

  area <- 0
  for(i in 1:np) {
    ii <- c(i, nextp[i])
    xpii <- xp[ii]
    ypii <- yp[ii]
    for(j in 1:nq) {
      jj <- c(j, nextq[j])
      area <- area +
        overlap.trapezium(xpii, ypii, xq[jj], yq[jj])
    }
  }
  return(area)
}

overlap.trapezium <- function(xa, ya, xb, yb, verb=F) {
  # compute area of overlap of two trapezia
  # which have same baseline y = 0
  #
  # first trapezium has vertices
  # (xa[1], 0), (xa[1], ya[1]), (xa[2], ya[2]), (xa[2], 0).
  # Similarly for second trapezium
  
  # Test for vertical edges
  dxa <- diff(xa)
  dxb <- diff(xb)
  if(dxa == 0 || dxb == 0)
    return(0)

  # Order x coordinates, x0 < x1
  if(dxa > 0) {
    signa <- 1
    lefta <- 1
    righta <- 2
    if(verb) cat("A is positive\n")
  } else {
    signa <- -1
    lefta <- 2
    righta <- 1
    if(verb) cat("A is negative\n")
  }
  if(dxb > 0) {
    signb <- 1
    leftb <- 1
    rightb <- 2
    if(verb) cat("B is positive\n")
  } else {
    signb <- -1
    leftb <- 2
    rightb <- 1
    if(verb) cat("B is negative\n")
  }
  signfactor <- signa * signb # actually (-signa) * (-signb)
  if(verb) cat(paste("sign factor =", signfactor, "\n"))

  # Intersect x ranges
  x0 <- max(xa[lefta], xb[leftb])
  x1 <- min(xa[righta], xb[rightb])
  if(x0 >= x1)
    return(0)
  if(verb) {
    cat(paste("Intersection of x ranges: [", x0, ",", x1, "]\n"))
    abline(v=x0, lty=3)
    abline(v=x1, lty=3)
  }

  # Compute associated y coordinates
  slopea <- diff(ya)/diff(xa)
  y0a <- ya[lefta] + slopea * (x0-xa[lefta])
  y1a <- ya[lefta] + slopea * (x1-xa[lefta])
  slopeb <- diff(yb)/diff(xb)
  y0b <- yb[leftb] + slopeb * (x0-xb[leftb])
  y1b <- yb[leftb] + slopeb * (x1-xb[leftb])
  
  # Determine whether upper edges intersect
  # if not, intersection is a single trapezium
  # if so, intersection is a union of two trapezia

  yd0 <- y0b - y0a
  yd1 <- y1b - y1a
  if(yd0 * yd1 >= 0) {
    # edges do not intersect
    areaT <- (x1 - x0) * (min(y1a,y1b) + min(y0a,y0b))/2
    if(verb) cat(paste("Edges do not intersect\n"))
  } else {
    # edges do intersect
    # find intersection
    xint <- x0 + (x1-x0) * abs(yd0/(yd1 - yd0))
    yint <- y0a + slopea * (xint - x0)
    if(verb) {
      cat(paste("Edges intersect at (", xint, ",", yint, ")\n"))
      points(xint, yint, cex=2, pch="O")
    }
    # evaluate left trapezium
    left <- (xint - x0) * (min(y0a, y0b) + yint)/2
    # evaluate right trapezium
    right <- (x1 - xint) * (min(y1a, y1b) + yint)/2
    areaT <- left + right
    if(verb)
      cat(paste("Left area = ", left, ", right=", right, "\n"))    
  }

  # return area of intersection multiplied by signs 
  return(signfactor * areaT)
}

