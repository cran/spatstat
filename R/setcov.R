#
#
#     setcov.R
#
#     $Revision: 1.2 $ $Date: 2002/07/18 10:54:34 $
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
  # Now set up a raster image 
  xstep <- mW$xstep
  ystep <- mW$ystep
  out <- im(G, xcol=xstep * ((-nc):nc), yrow=ystep * ((-nr):nr))
  return(out)
}

