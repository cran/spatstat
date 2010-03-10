#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 5.31 $	$Date: 2010/03/08 08:23:04 $
#
#
# -------- functions ----------------------------------------
#	Kcross()	cross-type K function K_{ij}
#                       between types i and j
#
#	Kdot()          K_{i\bullet}
#                       between type i and all points regardless of type
#
#       Kmulti()        (generic)
#
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#				including 'marks' vector
#	r		distance values at which to compute K	
#
# -------- standard output ------------------------------
#      A data frame with columns named
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

"Lcross" <- function(X, i, j, ...) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- levels(marks(X))[1]
  if(missing(j)) j <- levels(marks(X))[2]
  K <- Kcross(X, i, j, ...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  L <- rebadge.fv(L,
                  substitute(L[i,j](r),
                             list(i=paste(i),j=paste(j))),
                  "Lcross",
                  new.yexp=substitute(L[list(i,j)](r),
                                      list(i=paste(i),j=paste(j))))
  return(L)  
}

"Ldot" <- function(X, i, ...) {
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(i)) i <- levels(marks(X))[1]
  K <- Kdot(X, i, ...)
  L <- eval.fv(sqrt(K/pi))
  # relabel the fv object
  L <- rebadge.fv(L,
                  substitute(Ldot[i](r), list(i=paste(i))),
                  "Ldot")
  return(L)  
}

"Kcross" <- 
function(X, i, j, r=NULL, breaks=NULL,
         correction =c("border", "isotropic", "Ripley", "translate") , ...)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL
  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
  if(missing(j))
    j <- levels(marx)[2]
  I <- (marx == i)
  if(!any(I))
    stop(paste("No points have mark i =", i))

  if(i == j) {
    result <- Kest(X[I],
                   r=r, breaks=breaks, correction=correction, ...)
  } else {
    J <- (marx == j)
    if(!any(J))
      stop(paste("No points have mark j =", j))
    result <- Kmulti(X, I, J,
                     r=r, breaks=breaks, correction=correction, ...)
  }
  result <-
    rebadge.fv(result, 
               substitute(Kcross[i,j](r), list(i=paste(i),j=paste(j))),
               "Kcross",
               new.yexp=substitute(K[list(i,j)](r),
                                   list(i=paste(i),j=paste(j))))
  return(result)
}

"Kdot" <- 
function(X, i, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...)
{
  verifyclass(X, "ppp")
  if(!is.multitype(X, dfok=FALSE)) 
	stop("Point pattern must be multitype")
  if(missing(correction))
    correction <- NULL

  marx <- marks(X)
  if(missing(i))
    i <- levels(marx)[1]
        
  I <- (marx == i)
  J <- rep(TRUE, X$n)  # i.e. all points
	
  if(!any(I)) stop(paste("No points have mark i =", i))
	
  result <- Kmulti(X, I, J,
                   r=r, breaks=breaks, correction=correction, ...)
  result <-
    rebadge.fv(result,
               substitute(Kdot[i](r), list(i=paste(i))),
               "Kdot")
  return(result)
}


"Kmulti"<-
function(X, I, J, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...)
{
  verifyclass(X, "ppp")

  npoints <- X$n
  x <- X$x
  y <- X$y
  W <- X$window
  area <- area.owin(W)

  correction.given <- !missing(correction) && !is.null(correction)
  if(is.null(correction))
    correction <- c("border", "isotropic", "Ripley", "translate")
  
  correction <- pickoption("correction", correction,
                           c(none="none",
                             border="border",
                             "bord.modif"="bord.modif",
                             isotropic="isotropic",
                             Ripley="isotropic",
                             translate="translate",
                             best="best"),
                           multi=TRUE)

  correction <- implemented.for.K(correction, W$type, correction.given)

  if(!is.logical(I) || !is.logical(J))
    stop("I and J must be logical vectors")
  if(length(I) != npoints || length(J) != npoints)
    stop(paste("The length of I and J must equal",
               "the number of points in the pattern"))
	
  if(!any(I)) stop("no points satisfy I")
  if(!any(J)) stop("no points satisfy J")
		
  nI <- sum(I)
  nJ <- sum(J)
  lambdaI <- nI/area
  lambdaJ <- nJ/area

  # r values 
  rmaxdefault <- rmax.rule("K", W, lambdaJ)
  breaks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
        
  # recommended range of r values
  alim <- c(0, min(rmax, rmaxdefault))
        
  # this will be the output data frame
  # It will be given more columns later
  K <- data.frame(r=r, theo= pi * r^2)
  desc <- c("distance argument r", "theoretical Poisson %s")
  K <- fv(K, "r", substitute(Kmulti(r), NULL),
          "theo", , alim, c("r","%s[pois](r)"), desc, fname="Kmulti")

  # find close pairs of points
  XI <- X[I]
  XJ <- X[J]
  close <- crosspairs(XI, XJ, max(r))
# close$i and close$j are serial numbers in XI and XJ respectively;        
# map them to original serial numbers in X
  orig <- seq(npoints)
  imap <- orig[I]
  jmap <- orig[J]
  iX <- imap[close$i]
  jX <- jmap[close$j]
# eliminate any identical pairs
  if(any(I & J)) {
    ok <- (iX != jX)
    if(!all(ok)) {
      close$i  <- close$i[ok]
      close$j  <- close$j[ok]
      close$xi <- close$xi[ok]
      close$yi <- close$yi[ok]
      close$xj <- close$xj[ok]
      close$yj <- close$yj[ok]
      close$dx <- close$dx[ok]
      close$dy <- close$dy[ok]
      close$d  <- close$d[ok]
    }
  }
# extract information for these pairs (relative to orderings of XI, XJ)
  dcloseIJ <- close$d
  icloseI  <- close$i
  jcloseJ  <- close$j
        
# Compute estimates by each of the selected edge corrections.
        
  if(any(correction == "none")) {
    # uncorrected! 
    wh <- whist(dcloseIJ, breaks$val)  # no weights
    Kun <- cumsum(wh)/(lambdaI * lambdaJ * area)
    K <- bind.fv(K, data.frame(un=Kun), "%s[un](r)",
                 "uncorrected estimate of %s",
                 "un")
  }
  if(any(correction == "border" | correction == "bord.modif")) {
    # border method
    # distance to boundary from each point of type I
    bI <- bdist.points(XI)
    # distance to boundary from first element of each (i, j) pair
    bcloseI <- bI[icloseI]
    # apply reduced sample algorithm
    RS <- Kount(dcloseIJ, bcloseI, bI, breaks)
    if(any(correction == "bord.modif")) {
      denom.area <- eroded.areas(W, r)
      Kbm <- RS$numerator/(denom.area * nI * nJ)
      K <- bind.fv(K, data.frame(bord.modif=Kbm), "%s[bordm](r)",
                   "modified border-corrected estimate of %s",
                   "bord.modif")
    }
    if(any(correction == "border")) {
      Kb <- RS$numerator/(lambdaJ * RS$denom.count)
      K <- bind.fv(K, data.frame(border=Kb), "%s[bord](r)",
                   "border-corrected estimate of %s",
                   "border")
    }
  }
  if(any(correction == "translate")) {
    # translation correction
    edgewt <- edge.Trans(XI[icloseI], XJ[jcloseJ], paired=TRUE)
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    Ktrans <- cumsum(wh)/(lambdaI * lambdaJ * area)
    rmax <- diameter(W)/2
    Ktrans[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(trans=Ktrans), "%s[trans](r)",
                 "translation-corrected estimate of %s",
                 "trans")
  }
  if(any(correction == "isotropic")) {
    # Ripley isotropic correction
    edgewt <- edge.Ripley(XI[icloseI], matrix(dcloseIJ, ncol=1))
    wh <- whist(dcloseIJ, breaks$val, edgewt)
    Kiso <- cumsum(wh)/(lambdaI * lambdaJ * area)
    rmax <- diameter(W)/2
    Kiso[r >= rmax] <- NA
    K <- bind.fv(K, data.frame(iso=Kiso), "%s[iso](r)",
                 "Ripley isotropic correction estimate of %s",
                 "iso")
  }
  # default is to display them all
  attr(K, "fmla") <- . ~ r
  unitname(K) <- unitname(X)
  return(K)
}
