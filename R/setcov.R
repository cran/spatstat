#
#
#     setcov.R
#
#     $Revision: 1.1 $ $Date: 2002/05/06 09:39:58 $
#
#    Compute the set covariance function of a window
#
#

setcov <- function(W) {
  W <- as.owin(W)
  # pixel approximation
  mW <- as.mask(W)
  M <- mW$m
  # pad with zeroes
  nr <- nrow(M)
  nc <- ncol(M)
  Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Mpad[1:nr, 1:nc] <- M
  lengthMpad <- 4 * nc * nr
  # compute set covariance by fft
  fM <- fft(Mpad)
  G <- fft(Mod(fM)^2, inverse=TRUE)/lengthMpad
#  cat(paste("maximum imaginary part=", max(Im(G)), "\n"))
  G <- Mod(G) * mW$xstep * mW$ystep
  # Currently G[i,j] corresponds to a vector shift of
  #     dy = (i-1) mod nr, dx = (j-1) mod nc.
  # Rearrange this periodic function so that 
  # the origin of translations (0,0) is at matrix position (nr,nc)
  G <- G[ ((-nr):nr) %% (2 * nr) + 1, (-nc):nc %% (2*nc) + 1]
  # Now set up a raster image structure
  # using conventions similar to those for "owin" windows of type "mask"
  xstep <- mW$xstep
  ystep <- mW$ystep
  out <- list(m=G,
              dim=dim(G),
              xstep=xstep,
              ystep=ystep,
              xcol=xstep * ((-nc):nc),
              yrow=ystep * ((-nr):nr),
              xrange=xstep * nc * c(-1,1),
              yrange=ystep * nr * c(-1,1))
  class(out) <- "im"
  return(out)
}

# ----------------------------------------------------------
# Temporary code until we sort out the class structure
#
#
# This function is similar to nearest.raster.point except for
# the third argument 'im' and the different idiom for calculating
# row & column - which could be used in nearest.raster.point()

nearest.pixel <- function(x,y,im) {
  verifyclass(im, "im")
  nr <- im$dim[1]
  nc <- im$dim[2]
  cc <- round(1 + (x - im$xcol[1])/im$xstep)
  rr <- round(1 + (y - im$yrow[1])/im$ystep)
  cc <- pmax(1,pmin(cc, nc))
  rr <- pmax(1,pmin(rr, nr))
  return(list(row=rr, col=cc))
}

# This function is a generalisation of inside.owin()
# to images other than binary-valued images.

lookup.im <- function(im, x, y) {
  verifyclass(im, "im")

  if(length(x) != length(y))
    stop("x and y must be numeric vectors of equal length")
  value <- rep(NA, length(x))
               
  # test whether inside bounding rectangle
  xr <- range(im$xcol)
  yr <- range(im$yrow)
  frameok <- (xr[1] <= x) & (x <= xr[2]) & (yr[1] <= y) & (y <= yr[2])
  value[!frameok] <- 0
  
  if(all(!frameok))  # all points OUTSIDE range - no further work needed
    return(value)  # all zero

  # consider only those points which are inside the frame
  xf <- x[frameok]
  yf <- y[frameok]
  # map locations to raster (row,col) coordinates
  loc <- nearest.pixel(xf,yf,im)
  # look up mask values
  m <- im$m
  nf <- sum(frameok)
  vf <- logical(nf)
  for(i in 1:nf) 
    vf[i] <- m[loc$row[i],loc$col[i]]
  
  # insert into 'ok' vector
  value[frameok] <- vf

  if(any(is.na(value)))
    warning("Internal error: NA's generated")
  
  return(value)
}
  
