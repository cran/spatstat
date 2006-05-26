#
#	Kest.S		Estimation of K function
#
#	$Revision: 5.22 $	$Date: 2006/05/25 09:59:05 $
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

"Kest"<-
function(X, ..., r=NULL, breaks=NULL, 
         correction=c("border", "isotropic", "Ripley", "translate"),
         nlarge=500)
{
  verifyclass(X, "ppp")
  nlarge.given <- !missing(nlarge)
  rfixed <- !missing(r) || !missing(breaks)
        
  npoints <- X$n
  W <- X$window
  area <- area.owin(W)
  lambda <- npoints/area
  lambda2 <- (npoints * (npoints - 1))/(area^2)

  breaks <- handle.r.b.args(r, breaks, W)
  r <- breaks$r

  correction.given <- !missing(correction)
  correction.name <- c("border", "bord.modif", "isotropic", "Ripley", "translate")
  correction.list <- c("border", "bord.modif", "isotropic", "isotropic", "translate")
  correction.id <- pmatch(correction, correction.name)
  if(any(nbg <- is.na(correction.id)))
    stop(paste("unrecognised correction",
               if(sum(nbg) > 1) "s",
               ": ",
               paste(correction[nbg], collapse=", "),
               sep=""))
  correction <- correction.list[correction.id]

  # available selection of edge corrections depends on window
  if(W$type != "rectangle") {
    iso <- (correction == "isotropic") 
    if(any(iso)) {
      if(correction.given)
        warning("Isotropic correction not implemented for non-rectangular windows")
      correction <- correction[!iso]
    }
  }

  # recommended range of r values
  alim <- c(0, min(max(r),
                   diff(W$xrange)/4,
                   diff(W$yrange)/4,
                   sqrt(1000/(lambda * pi))))

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
    if(!rfixed) {
      rmax <- min(max(r), alim[2])
      r <- seq(0, rmax, length=length(r))
    }
    return(Kborder.engine(X, max(r), length(r), correction))
  }

  ###########################################
  # Slower code
  ###########################################
        
        
  # this will be the output data frame
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson K(r)")
  K <- fv(K, "r", substitute(K(r), NULL),
          "theo", , alim, c("r","Kpois(r)"), desc)

  # matrix of pairwise distances
  d <- pairdist(X$x, X$y)

  if(any(correction == "border" | correction == "bord.modif")) {
  # border method
  # Compute distances to boundary
    b <- bdist.points(X)
  # Ignore pairs (i,i)
    diag(d) <- Inf
  # apply reduced sample algorithm
    RS <- Kount(d, b, breaks, slow=FALSE)
    if(any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      Kbm <- RS$numerator/(lambda2 * denom.area)
      K <- bind.fv(K, data.frame(bord.modif=Kbm), "Kbord*(r)",
                   "modified border-corrected estimate of K(r)",
                   "bord.modif")
    }
    if(any(correction == "border")) {
      Kb <- RS$numerator/(lambda * RS$denom.count)
      K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                   "border-corrected estimate of K(r)",
                   "border")
    }
    # reset diagonal to original values
    diag(d) <- 0
  }

  if(any(c("translate", "isotropic") %in% correction)) {
    offdiag <- (row(d) != col(d))
    d.offdiag <- d[offdiag]
  }
  
  if(any(correction == "translate")) {
    # translation correction
    edgewt <- edge.Trans(X)
    wh <- whist(d.offdiag, breaks$val, edgewt[offdiag])
    Ktrans <- cumsum(wh)/(lambda2 * area)
    rmax <- diameter(W)/2
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                 "translation-corrected estimate of K(r)",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(X, d)
    wh <- whist(d.offdiag, breaks$val, edgewt[offdiag])
    Kiso <- cumsum(wh)/(lambda2 * area)
    rmax <- diameter(W)/2
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "Kiso(r)",
                 "Ripley isotropic correction estimate of K(r)",
                 "iso")
  }
   

  # default is to display them all
  attr(K, "fmla") <- . ~ r
  return(K)
}

################################################################  
#############  SUPPORTING ALGORITHMS ###########################
################################################################  

Kount <- function(d, b, breaks, slow=FALSE) {
  #
  # "internal" routine to compute border-correction estimate of K or Kij
  #
  # d : matrix of pairwise distances
  #                  (to exclude diagonal entries, set diag(d) = Inf)
  # b : column vector of distances to window boundary
  # breaks : breakpts object
  #

  if(slow) { ########## slow ##############
          
       r <- breaks$r
       
       nr <- length(r)
       numerator <- numeric(nr)
       denom.count <- numeric(nr)

       for(i in 1:nr) {
         close <- (d <= r[i])
         nclose <- matrowsum(close) # assumes diag(d) set to Inf
         bok <- (b > r[i])
         numerator[i] <- sum(nclose[bok])
         denom.count[i] <- sum(bok)
       }
	
  } else { ############# fast ####################

        # determine which distances d_{ij} were observed without censoring
        bb <- matrix(b, nrow=nrow(d), ncol=ncol(d))
        uncen <- (d <= bb)
        #
        # histogram of noncensored distances
        nco <- whist(d[uncen], breaks$val)
        # histogram of censoring times for noncensored distances
        ncc <- whist(bb[uncen], breaks$val)
        # histogram of censoring times (yes, this is a different total size)
        cen <- whist(b, breaks$val)
        # go
        RS <- reduced.sample(nco, cen, ncc, show=TRUE)
        # extract results
        numerator <- RS$numerator
        denom.count <- RS$denominator
        # check
        if(length(numerator) != breaks$ncells)
          stop("internal error: length(numerator) != breaks$ncells")
        if(length(denom.count) != breaks$ncells)
          stop("internal error: length(denom.count) != breaks$ncells")
  }
  
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
  desc <- c("distance argument r", "theoretical Poisson K(r)")
  Kfv <- fv(Kdf, "r", substitute(K(r), NULL),
          "theo", , c(0,rmax), c("r","Kpois(r)"), desc)

  ####### start computing ############
  # sort in ascending order of x coordinate
  Xsort <- X[order(X$x)]
  x <- Xsort$x
  y <- Xsort$y
  
  # boundary distances
  b <- bdist.points(Xsort)

  # call the C code
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
              PACKAGE="spatstat")
    area <- area.owin(W)
    if("bord.modif" %in% correction) {
      lambda2 <- (npoints * (npoints - 1))/(area^2)
      denom.area <- eroded.areas(W, r)
      bordmod <- res$numer/(lambda2 * denom.area)
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bordmod), "Kbord*(r)",
                   "modified border-corrected estimate of K(r)",
                   "bord.modif")
    }
    if("border" %in% correction) {
      lambda <- npoints/area
      bord <- res$numer/(lambda * res$denom)
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "Kbord(r)",
                   "border-corrected estimate of K(r)",
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
    res <- .C("Kwborder",
              nxy=as.integer(npoints),
              x=as.double(x),
              y=as.double(y),
              w=as.double(weights),
              b=as.double(b),
              nr=as.integer(nr),
              rmax=as.double(rmax),
              numer=as.double(numeric(nr)),
              denom=as.double(numeric(nr)),
              PACKAGE="spatstat")
    if("border" %in% correction) {
      bord <- res$numer/res$denom
      Kfv <- bind.fv(Kfv, data.frame(border=bord), "Kbord(r)",
                     "border-corrected estimate of Kinhom(r)",
                     "border")
    }
    if("bord.modif" %in% correction) {
      bm <- res$numer/eroded.areas(W, r)
      Kfv <- bind.fv(Kfv, data.frame(bord.modif=bm), "Kbord*(r)",
                     "modified border-corrected estimate of Kinhom(r)",
                     "bord.modif")
    }
  }
  ##
  # default is to display them all
  attr(Kfv, "fmla") <- . ~ r
  return(Kfv)
}

     

