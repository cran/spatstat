#
#	Kinhom.S	Estimation of K function for inhomogeneous patterns
#
#	$Revision: 1.8 $	$Date: 2005/06/08 18:09:25 $
#
#	Kinhom()	compute estimate of K_inhom
#
#       Currently uses border method and slow code...
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
#       lambda          either a vector of intensity values for points of X
#                       or a matrix of two-point conditional intensities
#                       for pairs of points of X
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
         correction=c("border", "isotropic", "Ripley", "translate"), ...)
{
    verifyclass(X, "ppp")
    W <- X$window
    npoints <- X$n
    area <- area.owin(W)
    breaks <- handle.r.b.args(r, breaks, X$window)
    r <- breaks$r
    if(is.vector(lambda)) {
      if (length(lambda) != npoints) 
        stop("The length of the vector \'lambda\' should equal the number of data points")
      weight <- 1/outer(lambda, lambda, "*")
    } else if(is.matrix(lambda)) {
      if(any(dim(lambda) != npoints))
        stop("The matrix \'lambda\' should be n x n where n is the number of data points")
      weight <- 1/lambda
    }

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
    K <- fv(K, "r", "K(r)", "theo", , alim, c("r","Kpois(r)"), desc)
        
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
      RS <- Kwtsum(d, b, weight, breaks, slow)
      if(any(correction == "border")) {
        Kb <- area * RS$numerator/RS$denom.sum
        K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                     "border-corrected estimate of K(r)",
                     "border")
      }
      if(any(correction == "bord.modif")) {
        denom.area <- eroded.areas(W, r)
        Kbm <- RS$numerator/denom.area
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


Kwtsum <- function(d, b, weight, breaks, slow=FALSE) {
  #
  # "internal" routine to compute border-correction estimate of Kinhom
  #
  # d : matrix of pairwise distances
  #                  (to exclude diagonal entries, set diag(d) = Inf)
  # b : column vector of distances to window boundary
  # weight: matrix of weights for x[i], x[j] pairs
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
  
  weightmargin <- matrowsum(weight)

  if(slow) { ########## slow ##############
          
       r <- breaks$r
       
       nr <- length(r)
       numerator <- numeric(nr)
       denom.sum <- numeric(nr)

       for(i in 1:nr) {
         close <- (d <= r[i])
         numer <- matrowsum(close * weight) # assumes diag(d) set to Inf
         bok <- (b > r[i])
         numerator[i] <- sum(numer[bok])
         denom.sum[i] <- sum(weightmargin[bok])
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
        numerator <- RS$numerator
        denom.sum <- RS$denominator
        # check
        if(length(numerator) != breaks$ncells)
          stop("internal error: length(numerator) != breaks$ncells")
        if(length(denom.sum) != breaks$ncells)
          stop("internal error: length(denom.count) != breaks$ncells")
  }
  
  return(list(numerator=numerator, denom.sum=denom.sum))
}
	
