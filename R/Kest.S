#
#	Kest.S		Estimation of K function
#
#	$Revision: 5.48 $	$Date: 2009/04/29 17:08:19 $
#
#
# -------- functions ----------------------------------------
#	Kest()		compute estimate of K
#                       using various edge corrections
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
# -------- standard output ------------------------------
#      A data frame (class "fv") with columns named
#
#	r:		same as input
#
#	trans:		K function estimated by translation correction
#
#	iso:		K function estimated by Ripley isotropic correction
#
#	theo:		K function for Poisson ( = pi * r ^2 )
#
#	border:		K function estimated by border method
#			using standard formula (denominator = count of points)
#
#       bord.modif:	K function estimated by border method
#			using modified formula 
#			(denominator = area of eroded window
#
# ------------------------------------------------------------------------

"Lest" <- function(...) {
  K <- Kest(...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  L <- rebadge.fv(L, substitute(L(r), NULL), "L")
  return(L)  
}

"Kest"<-
function(X, ..., r=NULL, breaks=NULL, 
         correction=c("border", "isotropic", "Ripley", "translate"),
         nlarge=3000, domain=NULL)
{
  verifyclass(X, "ppp")
  nlarge.given <- !missing(nlarge) && !is.null(nlarge)
  rfixed <- !is.null(r) || !is.null(breaks)
        
  npoints <- X$n
  W <- X$window
  area <- area.owin(W)
  lambda <- npoints/area
  lambda2 <- (npoints * (npoints - 1))/(area^2)

  if(!is.null(domain)) {
    # estimate based on contributions from a subdomain
    domain <- as.owin(domain)
    if(!is.subset.owin(domain, X$window))
      stop(paste(dQuote("domain"),
                 "is not a subset of the window of X"))
    # trick Kdot() into doing it
    indom <- inside.owin(X$x, X$y, domain)
    Kd <- Kdot(X %mark% indom, i="TRUE",
               r=r, breaks=breaks, correction=correction)
    # correct for estimate of intensity of type i
    adj <- area/area.owin(domain)
    Kd <- eval.fv(adj * Kd)
    # relabel and exit
    Kd <- rebadge.fv(Kd, substitute(K(r), NULL), "K")
    return(Kd)
  }

  rmaxdefault <- rmax.rule("K", W, lambda)        
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  
  correction.given <- !missing(correction)
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate"),
                           multi=TRUE)

  # available selection of edge corrections depends on window
  if(W$type == "mask") {
    iso <- (correction == "isotropic") 
    if(any(iso)) {
      if(correction.given)
        warning("Isotropic correction not implemented for binary masks")
      correction <- correction[!iso]
    }
  }

  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))

  ###########################################
  # Efficient code for border method
  # Usable only if r values are evenly spaced from 0 to rmax
  # Invoked automatically if number of points is large

  can.do.fast <- breaks$even
  large.n    <- (npoints >= nlarge)
  borderonly <- all(correction == "border" | correction == "bord.modif")
  will.do.fast <- can.do.fast && (borderonly || large.n)
  asked      <- borderonly || (nlarge.given && large.n)
  if(will.do.fast && !asked)
        message(paste("number of data points exceeds",
                      nlarge, "- computing border estimate only"))
  if(asked && !can.do.fast)
    warning("r values not evenly spaced - cannot use efficient code")

  if(will.do.fast) {
    # restrict r values to recommended range, unless specifically requested
    if(!rfixed) 
      r <- seq(0, alim[2], length=length(r))
    return(Kborder.engine(X, max(r), length(r), correction))
  }

  ###########################################
  # Slower code
  ###########################################
        
        
  # this will be the output data frame
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", substitute(K(r), NULL),
          "theo", , alim, c("r","%spois(r)"), desc, fname="K")

  # identify all close pairs
  rmax <- max(r)
  close <- closepairs(X, rmax)
  DIJ <- close$d
  XI <- ppp(close$xi, close$yi, window=W, check=FALSE)
  
  if(any(correction == "none")) {
    # uncorrected! For demonstration purposes only!
    wh <- whist(DIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/(lambda2 * area)
    K <- bind.fv(K, data.frame(un=Kun), "%sun(r)",
                 "uncorrected estimate of %s",
                 "un")
  }
  
  if(any(correction == "border" | correction == "bord.modif")) {
  # border method
  # Compute distances to boundary
    b <- bdist.points(X)
    I <- close$i
    bI <- b[I]
  # apply reduced sample algorithm
    RS <- Kount(DIJ, bI, b, breaks)
    if(any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      Kbm <- RS$numerator/(lambda2 * denom.area)
      K <- bind.fv(K, data.frame(bord.modif=Kbm), "%sbord*(r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
    }
    if(any(correction == "border")) {
      Kb <- RS$numerator/(lambda * RS$denom.count)
      K <- bind.fv(K, data.frame(border=Kb), "%sbord(r)",
                   "border-corrected estimate of %s",
                   "border")
    }
  }

  if(any(correction == "translate")) {
    # translation correction
    XJ <- ppp(close$xj, close$yj, window=W, check=FALSE)
    edgewt <- edge.Trans(XI, XJ, paired=TRUE)
    wh <- whist(DIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/(lambda2 * area)
    h <- diameter(W)/2
    Ktrans[r >= h] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "%strans(r)",
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI, matrix(DIJ, ncol=1))
    wh <- whist(DIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/(lambda2 * area)
    h <- diameter(W)/2
    Kiso[r >= h] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "%siso(r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
  }
   

  # default is to display them all
  attr(K, "fmla") <- . ~ r
  unitname(K) <- unitname(X)
  return(K)
}

################################################################  
#############  SUPPORTING ALGORITHMS ###########################
################################################################  

Kount <- function(dIJ, bI, b, breaks) {
  #
  # "internal" routine to compute border-correction estimate of K or Kij
  #
  # dIJ:  vector containing pairwise distances for selected I,J pairs
  # bI:   corresponding vector of boundary distances for I
  # b:    vector of ALL distances to window boundary
  #
  # breaks : breakpts object
  #

  stopifnot(length(dIJ) == length(bI))
  
  # determine which distances d_{ij} were observed without censoring
  uncen <- (dIJ <= bI)
  # histogram of noncensored distances
  nco <- whist(dIJ[uncen], breaks$val)
  # histogram of censoring times for noncensored distances
  ncc <- whist(bI[uncen], breaks$val)
  # histogram of censoring times (yes, this is a different total size)
  cen <- whist(b, breaks$val)
  # count censoring times beyond rightmost breakpoint
  uppercen <- sum(b > max(breaks$val))
  # go
  RS <- reduced.sample(nco, cen, ncc, show=TRUE, uppercen=uppercen)
  # extract results
  numerator <- RS$numerator
  denom.count <- RS$denominator
  # check
  if(length(numerator) != breaks$ncells)
    stop("internal error: length(numerator) != breaks$ncells")
  if(length(denom.count) != breaks$ncells)
    stop("internal error: length(denom.count) != breaks$ncells")
  
  return(list(numerator=numerator, denom.count=denom.count))
}

#### interface to C code for border method

Kborder.engine <- function(X, rmax, nr=100,
                           correction=c("border", "bord.modif"),
                           weights) 
{
  verifyclass(X, "ppp")
  npoints <- X$n
  W <- X$window

  if(missing(rmax))
    rmax <- diameter(W)/4
  r <- seq(0, rmax, length=nr)

  # this will be the output data frame
  Kdf <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  Kfv <- fv(Kdf, "r", substitute(K(r), NULL),
          "theo", , c(0,rmax), c("r","%spois(r)"), desc, fname="K")

  ####### start computing ############
  # sort in ascending order of x coordinate
  orderX <- order(X$x)
  Xsort <- X[orderX]
  x <- Xsort$x
  y <- Xsort$y
  
  # boundary distances
  b <- bdist.points(Xsort)

  # call the C code
  DUP <- spatstat.options("dupC")
  
  if(missing(weights)) {
    res <- .C("Kborder",
              nxy=as.integer(npoints),
              x=as.double(x),
              y=as.double(y),
              b=as.double(b),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.integer(integer(nr)),
              denom=as.integer(integer(nr)),
              DUP=DUP,
              PACKAGE="spatstat")
    area <- area.owin(W)
    if("bord.modif" %in% correction) {
      lambda2 <- (npoints * (npoints - 1))/(area^2)
      denom.area <- eroded.areas(W, r)
      bordmod <- res$numer/(lambda2 * denom.area)
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bordmod), "%sbord*(r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
    }
    if("border" %in% correction) {
      lambda <- npoints/area
      bord <- res$numer/(lambda * res$denom)
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "%sbord(r)",
                   "border-corrected estimate of %s",
                   "border")
    }
  } else {
    # weighted version
    if(is.numeric(weights)) {
      if(length(weights) != X$n)
        stop("length of weights argument does not match number of points in X")
    } else {
      wim <- as.im(weights, W)
      weights <- wim[X, drop=FALSE]
      if(any(is.na(weights)))
        stop("domain of weights image does not contain all points of X")
    }
    weights.Xsort <- weights[orderX]
    res <- .C("Kwborder",
              nxy=as.integer(npoints),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights.Xsort),
              b=as.double(b),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)),
              denom=as.double(numeric(nr)),
              DUP=DUP,
              PACKAGE="spatstat")
    if("border" %in% correction) {
      bord <- res$numer/res$denom
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "%sbord(r)",
                     "border-corrected estimate of %sinhom(r)",
                     "border")
    }
    if("bord.modif" %in% correction) {
      bm <- res$numer/eroded.areas(W, r)
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bm), "%sbord*(r)",
                     "modified border-corrected estimate of %s",
                     "bord.modif")
    }
  }
  ##
  # default is to display them all
  attr(Kfv, "fmla") <- . ~ r
  unitname(Kfv) <- unitname(X)
  return(Kfv)
}

     

rmax.rule <- function(fun="K", W, lambda) {
  verifyclass(W, "owin")
  switch(fun,
         K = {
           # Ripley's Rule
           ripley <- min(diff(W$xrange), diff(W$yrange))/4
           # Count at most 1000 neighbours per point
           rlarge <- if(!missing(lambda)) sqrt(1000 /(pi * lambda)) else Inf
           rmax <- min(rlarge, ripley)
         },
         F = ,
         G = ,
         J = {
           # rule of thumb
           rdiam  <- diameter(W)/2
           # Poisson process has F(rlarge) = 1 - 10^(-5)
           rlarge <-
             if(!missing(lambda)) sqrt(log(10^5)/(pi * lambda)) else Inf
           rmax <- min(rlarge, rdiam)
         },
         stop(paste("Unrecognised function type", sQuote(fun)))
         )
  return(rmax)
}
           
           
