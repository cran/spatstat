#
#	localK.R		Getis-Franklin neighbourhood density function
#
#	$Revision: 1.13 $	$Date: 2009/07/24 06:51:58 $
#
#

"localL" <-
  function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  localK(X, wantL=TRUE,
         correction=correction, verbose=verbose, rvalue=rvalue)
}

"localK" <-
  function(X, ..., correction="Ripley", verbose=TRUE, rvalue=NULL)
{
  verifyclass(X, "ppp")
  
  wantL <- resolve.defaults(list(...), list(wantL = FALSE))$wantL

  npoints <- X$n 
  W <- X$window
  area <- area.owin(W)
  lambda <- npoints/area
  lambda1 <- (npoints - 1)/area
  lambda2 <- (npoints * (npoints - 1))/(area^2)

  if(is.null(rvalue)) 
    rmaxdefault <- rmax.rule("K", W, lambda)
  else {
    stopifnot(is.numeric(rvalue))
    stopifnot(length(rvalue) == 1)
    stopifnot(rvalue >= 0)
    rmaxdefault <- rvalue
  }
  breaks <- handle.r.b.args(NULL, NULL, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  
  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             best="best"),
                           multi=FALSE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  # identify all close pairs
  rmax <- max(r)
  close <- closepairs(X, rmax)
  DIJ <- close$d
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
  I <- close$i

  # initialise
  df <- as.data.frame(matrix(NA, length(r), npoints))
  labl <- desc <- character(npoints)
  fname <- if(wantL) "L(r)" else "K(r)"

  bkt <- function(x) { paste("[", x, "]", sep="") }
  
  switch(correction,
         none={
           # uncorrected! For demonstration purposes only!
           ftag <- if(wantL) "L[un]" else "K[un]"
           for(i in 1:npoints) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val)  # no weights
             df[,i] <- cumsum(wh)/lambda1
             icode <- numalign(i, npoints)
             names(df)[i] <- paste("un", icode, sep="")
             labl[i] <- paste(ftag, bkt(icode), "(r)", sep="")
             desc[i] <- paste("uncorrected estimate of", fname,
                              "for point", icode)
             if(verbose) progressreport(i, npoints)
           }
         },
         translate={
           # Translation correction
           ftag <- if(wantL) "L[trans]" else "K[trans]"
           XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
           edgewt <- edge.Trans(XI, XJ, paired=TRUE)
           for(i in 1:npoints) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
             Ktrans <- cumsum(wh)/lambda1
             df[,i] <- Ktrans
             icode <- numalign(i, npoints)
             names(df)[i] <- paste("trans", icode, sep="")
             labl[i] <- paste(ftag, bkt(icode), "(r)", sep="")
             desc[i] <- paste("translation-corrected estimate of", fname,
                              "for point", icode)
             if(verbose) progressreport(i, npoints)
           }
           h <- diameter(W)/2
           df[r >= h, ] <- NA
         },
         isotropic={
           # Ripley isotropic correction
           ftag <- if(wantL) "L[iso]" else "K[iso]"
           edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
           for(i in 1:npoints) {
             ii <- (I == i)
             wh <- whist(DIJ[ii], breaks$val, edgewt[ii])
             Kiso <- cumsum(wh)/lambda1
             df[,i] <- Kiso
             icode <- numalign(i, npoints)
             names(df)[i] <- paste("iso", icode, sep="")
             labl[i] <- paste(ftag, bkt(icode), "(r)", sep="")
             desc[i] <- paste("Ripley isotropic correction estimate of", fname,
                              "for point", icode)
             if(verbose) progressreport(i, npoints)
           }
           h <- diameter(W)/2
           df[r >= h, ] <- NA
         })
  # transform values if L required
  if(wantL)
    df <- sqrt(df/pi)
  
  # return vector of values at r=rvalue, if desired
  if(!is.null(rvalue)) {
    nr <- length(r)
    if(r[nr] != rvalue)
      stop("Internal error - rvalue not attained")
    return(as.numeric(df[nr,]))
  }
  # function value table required
  # add r and theo
  if(!wantL) {
    ylab <- substitute(localK(r), NULL)
    df <- cbind(df, data.frame(r=r, theo=pi * r^2))
    fnam <- "localK"
  } else {
    ylab <- substitute(localL(r), NULL)
    df <- cbind(df, data.frame(r=r, theo=r))
    fnam <- "localL"
  }
  desc <- c(desc, c("distance argument r", "theoretical Poisson %s"))
  labl <- c(labl, c("r", "%s[pois](r)"))
  # create fv object
  K <- fv(df, "r", ylab, "theo", , alim, labl, desc, fname=fnam)
  # default is to display them all
  attr(K, "fmla") <- . ~ r
  unitname(K) <- unitname(X)
  return(K)
}


