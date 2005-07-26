#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.12 $	$Date: 2005/07/21 10:22:34 $
#
#	Kinhom()	compute estimate of K_inhom
#
#
#       Reference:
#            Non- and semiparametric estimation of interaction
#	     in inhomogeneous point patterns
#            A.Baddeley, J.Moller, R.Waagepetersen
#            Statistica Neerlandica 54 (2000) 329--350.
#
# -------- functions ----------------------------------------
#	Kinhom()	compute estimate of K
#                       using various edge corrections
#
#       Kwtsum()         internal routine for border correction
#
# -------- standard arguments ------------------------------	
#	X		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
#       lambda          vector of intensity values for points of X
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
#			(denominator = sum of weights of points)
#
#       bord.modif:	K function estimated by border method
#			(denominator = area of eroded window)
#
# ------------------------------------------------------------------------

"Kinhom"<-
  function (X, lambda, r = NULL, breaks = NULL, slow=FALSE,
         correction=c("border", "bord.modif", "isotropic", "translate"), ...,
            lambda2=NULL)
{
    verifyclass(X, "ppp")
    W <- X$window
    npoints <- X$n
    area <- area.owin(W)
    breaks <- handle.r.b.args(r, breaks, X$window)
    r <- breaks$r
    
    if(!is.vector(lambda))
      stop("\`lambda\' should be a vector")
    if(length(lambda) != npoints) 
      stop("The length of \'lambda\' should equal the number of data points")
    if(is.null(lambda2))
      weight <- 1/outer(lambda, lambda, "*")
    else {
      if(!is.matrix(lambda2))
        stop("lambda2 should be a matrix")
      if(nrow(lambda2) != ncol(lambda2))
        stop("lambda2 should be a square matrix")
      if(nrow(lambda2) != npoints)
        stop("Dimensions of lambda2 do not match number of data points")
      weight <- 1/lambda2
    } 
    
    # match corrections
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
      iso <- (correction == "isotropic") | (correction == "Ripley")
      if(any(iso)) {
        if(!missing(correction))
          warning("Isotropic correction not implemented for non-rectangular windows")
        correction <- correction[!iso]
      }
    }

    # recommended range of r values
    alim <- c(0, min(diff(X$window$xrange), diff(X$window$yrange))/4)
        
    # this will be the output data frame
    K <- data.frame(r=r, theo= pi * r^2)
    desc <- c("distance argument r", "theoretical Poisson K(r)")
    K <- fv(K, "r", substitute(Kinhom(r), NULL),
            "theo", , alim, c("r","Kpois(r)"), desc)
        
    # pairwise distance
    d <- pairdist(X$x, X$y)
    
    offdiag <- (row(d) != col(d))
        
    if(any(correction == "border" | correction == "bord.modif")) {
      # border method
      # Compute distances to boundary
      b <- bdist.points(X)
      # Ignore pairs (i,i)
      diag(d) <- Inf
      # apply reduced sample algorithm
      RS <- Kwtsum(d, b, weight, 1/lambda, breaks, slow)
      if(any(correction == "border")) {
        Kb <- RS$ratio
        K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                     "border-corrected estimate of K(r)",
                     "border")
      }
      if(any(correction == "bord.modif")) {
        Kbm <- RS$numerator/eroded.areas(W, r)
        K <- bind.fv(K, data.frame(bord.modif=Kbm), "Kbord*(r)",
                     "modified border-corrected estimate of K(r)",
                     "bord.modif")
      }
      # reset diagonal to original values
      diag(d) <- 0
    }
    if(any(correction == "translate")) {
      # translation correction
      edgewt <- edge.Trans(X)
      allweight <- edgewt * weight
      wh <- whist(d[offdiag], breaks$val, allweight[offdiag])
      Ktrans <- cumsum(wh)/area
      rmax <- diameter(W)/2
      Ktrans[r >= rmax] <- NA
      K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                   "translation-correction estimate of K(r)",
                   "trans")
    }
    if(any(correction == "isotropic" | correction == "Ripley")) {
      # Ripley isotropic correction
      edgewt <- edge.Ripley(X, d)
      allweight <- edgewt * weight
      wh <- whist(d[offdiag], breaks$val, allweight[offdiag])
      Kiso <- cumsum(wh)/area
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


Kwtsum <- function(d, b, weight, weightmargin, breaks, slow=FALSE) {
  #
  # "internal" routine to compute border-correction estimates of Kinhom
  #
  # d : matrix of pairwise distances
  #                  (to exclude diagonal entries, set diag(d) = Inf)
  # b : column vector of distances to window boundary
  # weight: matrix of weights for x[i], x[j] pairs
  # weightmargin: vector of weights for x[i]
  # breaks : breakpts object
  #

  if(!is.matrix(d))
    stop("\'d\' must be a matrix")
  if(!is.matrix(weight))
    stop("\'weight\' must be a matrix")
  if(any(dim(d) != dim(weight)))
    stop("matrices \'d\' and \'weight\' have different dimensions")
  if(length(b) != nrow(d))
    stop("length(b) does not match nrow(d)")
  
  if(slow) { ########## slow ##############
          
       r <- breaks$r
       
       nr <- length(r)
       numerator <- numeric(nr)
       denominator <- numeric(nr)

       for(i in 1:nr) {
         close <- (d <= r[i])
         numer <- matrowsum(close * weight) # assumes diag(d) set to Inf
         bok <- (b > r[i])
         numerator[i] <- sum(numer[bok])
         denominator[i] <- sum(weightmargin[bok])
       }
	
  } else { ############# fast ####################

        # determine which distances d_{ij} were observed without censoring
        bb <- matrix(b, nrow=nrow(d), ncol=ncol(d))
        uncen <- (d <= bb)
        #
        # histogram of noncensored distances
        nco <- whist(d[uncen], breaks$val, weight[uncen])
        # histogram of censoring times for noncensored distances
        ncc <- whist(bb[uncen], breaks$val, weight[uncen])
        # histogram of censoring times (yes, this is a different total size)
        cen <- whist(b, breaks$val, weightmargin)
        # go
        RS <- reduced.sample(nco, cen, ncc, show=TRUE)
        # extract results
        numerator   <- RS$numerator
        denominator <- RS$denominator
        ratio        <- RS$numerator/RS$denominator
        # check
        if(length(numerator) != breaks$ncells)
          stop("internal error: length(numerator) != breaks$ncells")
        if(length(denominator) != breaks$ncells)
          stop("internal error: length(denom.count) != breaks$ncells")
  }
  
  return(list(numerator=numerator, denominator=denominator, ratio=ratio))
}
	
