#
#
#     setcov.R
#
#     $Revision: 1.8 $ $Date: 2011/10/10 08:29:03 $
#
#    Compute the set covariance function of a window
#    or the (noncentred) spatial covariance function of an image
#

setcov <- function(W, V=W, ...) {
  W <- as.owin(W)
  # pixel approximation
  mW <- as.mask(W, ...)
  Z <- as.im(mW, na.replace=0)
  if(missing(V)) 
    return(imcov(Z))
  # cross-covariance
  V <- as.owin(V)
  mV <- as.mask(V, ...)
  Z2 <- as.im(mV, na.replace=0)
  imcov(Z, Z2)
}

imcov <- function(X, Y=X) {
  stopifnot(is.im(X))
  crosscov <- !missing(Y)
  if(crosscov) {
    # cross-covariance 
    stopifnot(is.im(Y))
    Xbox <- as.rectangle(X)
    Ybox <- as.rectangle(Y)
    # first shift images to common midpoint, to reduce storage
    Xmid <- centroid.owin(Xbox)
    Ymid <- centroid.owin(Ybox)
    svec <- as.numeric(Xmid) - as.numeric(Ymid)
    Y <- shift(Y, svec)
    # ensure images are compatible
    XY <- harmonise.im(X=X, Y=Y)
    X <- XY$X
    Y <- XY$Y
  }
  M <- X$v
  M[is.na(M)] <- 0
  xstep <- X$xstep
  ystep <- X$ystep
  # pad with zeroes
  nr <- nrow(M)
  nc <- ncol(M)
  Mpad <- matrix(0, ncol=2*nc, nrow=2*nr)
  Mpad[1:nr, 1:nc] <- M
  lengthMpad <- 4 * nc * nr
  fM <- fft(Mpad)
  if(!crosscov) {
    # compute set covariance by fft
    G <- fft(Mod(fM)^2, inverse=TRUE)/lengthMpad
  } else {
    # compute set cross-covariance by fft
    N <- Y$v
    N[is.na(N)] <- 0
    Npad <- matrix(0, ncol=2*nc, nrow=2*nr)
    Npad[1:nr, 1:nc] <- N
    fN <- fft(Npad)
    G <- fft(fM * Conj(fN), inverse=TRUE)/lengthMpad
  }
#  cat(paste("maximum imaginary part=", max(Im(G)), "\n"))
  G <- Mod(G) * xstep * ystep
  # Currently G[i,j] corresponds to a vector shift of
  #     dy = (i-1) mod nr, dx = (j-1) mod nc.
  # Rearrange this periodic function so that 
  # the origin of translations (0,0) is at matrix position (nr,nc)
  G <- G[ ((-nr):nr) %% (2 * nr) + 1, (-nc):nc %% (2*nc) + 1]
  # Now set up a raster image 
  out <- im(G, xcol=xstep * ((-nc):nc), yrow=ystep * ((-nr):nr))
  if(crosscov) {
    # restrict domain
    width <- diff(Xbox$xrange) + diff(Ybox$xrange)
    height <- diff(Xbox$yrange) + diff(Ybox$yrange)
    XYbox <- owin(c(-1,1) * width/2, c(-1,1) * height/2)
    out <- out[XYbox]
    # undo shift
    out <- shift(out, svec)
  }
  return(out)
}

