# place points at regular intervals along line segments

pointsOnLines <- function(X, eps=NULL, np=1000) {
  stopifnot(is.psp(X))
  win <- X$window
  len <- lengths.psp(X)
  nseg <- length(len)
  if(missing(eps)) {
    stopifnot(is.numeric(np) && length(np) == 1 & is.finite(np) && np > 0)
    eps <- sum(len)/np
  } else
  stopifnot(is.numeric(eps) && length(eps) == 1 && is.finite(eps) && eps > 0)
  Xdf    <- as.data.frame(X)
  Z <- ppp(numeric(0), numeric(0), window=win)
  for(i in 1:nseg) {
    # divide segment into pieces of length eps
    # with shorter bits at each end
    leni <- len[i]
    nwhole <- floor(leni/eps)
    if(leni/eps - nwhole < 0.5 && nwhole > 2)
      nwhole <- nwhole - 1
    rump <- (leni - nwhole * eps)/2
    brks <- c(0, rump + (0:nwhole) * eps, leni)
    nbrks <- length(brks)
    # points at middle of each piece
    ss <- (brks[-1] + brks[-nbrks])/2
    x <- with(Xdf, x0[i] + (ss/leni) * (x1[i]-x0[i]))
    y <- with(Xdf, y0[i] + (ss/leni) * (y1[i]-y0[i]))
    Zi <- ppp(x=x, y=y, window=win)
    Z <- superimpose(Z, Zi, W=win, check=FALSE)
  }
  return(Z)
}
