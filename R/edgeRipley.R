#
#        edgeRipley.R
#
#    $Revision: 1.1 $    $Date: 2002/04/07 11:14:43 $
#
#    Ripley isotropic edge correction weights
#
#  edge.Ripley(X, r, W)      compute isotropic correction weights
#                            for centres X[i], radii r[i,j], window W
#
#  To estimate the K-function see the idiom in "Kest.S"
#
#######################################################################

edge.Ripley <- function(X, r, W=X$window) {
  # X is a point pattern, or equivalent
  X <- as.ppp(X, W)
  W <- X$window
  if(W$type != "rectangle")
	stop("sorry, Ripley isotropic correction is only implemented\
for rectangular windows")

  if(!is.matrix(r) || nrow(r) != X$n)
    stop("r should be a matrix with nrow(r) = length(X$x)")

  x <- X$x
  y <- X$y

  # perpendicular distance from point to each edge of rectangle
  # L = left, R = right, D = down, U = up
  dL  <- x - W$xrange[1]
  dR  <- W$xrange[2] - x
  dD  <- y - W$yrange[1]
  dU  <- W$yrange[2] - y

  # detect whether any points are corners of the rectangle
  small <- function(x) { abs(x) < .Machine$double.eps }
  corner <- (small(dL) + small(dR) + small(dD) + small(dU) >= 2)
  
  # angle between (a) perpendicular to edge of rectangle
  # and (b) line from point to corner of rectangle
  bLU <- atan2(dU, dL)
  bLD <- atan2(dD, dL)
  bRU <- atan2(dU, dR)
  bRD <- atan2(dD, dR)
  bUL <- atan2(dL, dU)
  bUR <- atan2(dR, dU)
  bDL <- atan2(dL, dD)
  bDR <- atan2(dR, dD)

 # The above are all vectors [i]
 # Now we compute matrices [i,j]

  # half the angle subtended by the intersection between
  # the circle of radius r[i,j] centred on point i
  # and each edge of the rectangle (prolonged to an infinite line)

  hang <- function(d, r) {
    answer <- matrix(0, nrow=nrow(r), ncol=ncol(r))
    # replicate d[i] over j index
    d <- matrix(d, nrow=nrow(r), ncol=ncol(r))
    hit <- (d < r)
    answer[hit] <- acos(d[hit]/r[hit])
    answer
  }

  aL <- hang(dL, r)
  aR <- hang(dR, r)
  aD <- hang(dD, r) 
  aU <- hang(dU, r)

  # apply maxima
  # note: a* are matrices; b** are vectors;
  # b** are implicitly replicated over j index
  cL <- pmin(aL, bLU) + pmin(aL, bLD)
  cR <- pmin(aR, bRU) + pmin(aR, bRD)
  cU <- pmin(aU, bUL) + pmin(aU, bUR)
  cD <- pmin(aD, bDL) + pmin(aD, bDR)

  # total exterior angle
  ext <- cL + cR + cU + cD

  # add pi/2 for corners 
  if(any(corner))
    ext[corner,] <- ext[corner,] + pi/2

  # OK, now compute weight
  weight <- 1 / (1 - ext/(2 * pi))

  return(weight)
}
