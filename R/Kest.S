#
#	Kest.S		Estimation of K function
#
#	$Revision: 5.9 $	$Date: 2005/02/08 02:10:11 $
#
#
# -------- functions ----------------------------------------
#	Kest()		compute estimate of K
#                       using various edge corrections
#
#       Kount()         internal routine for border correction
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
function(X, r=NULL, breaks=NULL, slow=FALSE,
         correction=c("border", "isotropic", "Ripley", "translate"), ...)
{
	verifyclass(X, "ppp")

	npoints <- X$n
        W <- X$window
	area <- area.owin(W)
	lambda <- npoints/area
	lambda2 <- (npoints * (npoints - 1))/(area^2)

        breaks <- handle.r.b.args(r, breaks, W)
        r <- breaks$r

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
          RS <- Kount(d, b, breaks, slow)
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
        if(any(correction == "translate")) {
          # translation correction
            edgewt <- edge.Trans(X)
            wh <- whist(d[offdiag], breaks$val, edgewt[offdiag])
            Ktrans <- cumsum(wh)/(lambda2 * area)
            rmax <- diameter(W)/2
            Ktrans[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                         "translation-corrected estimate of K(r)",
                         "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(X, d)
            wh <- whist(d[offdiag], breaks$val, edgewt[offdiag])
            Kiso <- cumsum(wh)/(lambda2 * area)
            rmax <- diameter(W)/2
            Kiso[r >= rmax] <- NA
            K <- bind.fv(K, data.frame(iso=Kiso), "Kiso(r)",
                         "Ripley isotropic correction estimate of K(r)",
                         "iso")
        }

        # which corrections have been computed?
        nama2 <- names(K)
        corrxns <- rev(nama2[nama2 != "r"])

        # default is to display them all
        attr(K, "fmla") <- deparse(as.formula(paste(
                       "cbind(",
                        paste(corrxns, collapse=","),
                        ") ~ r")))
        return(K)
}
	
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
