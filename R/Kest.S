#
#	Kest.S		Estimation of K function
#
#	$Revision: 4.5 $	$Date: 2001/11/06 10:27:52 $
#
#	Border correction only.
#	Arbitrary window.
#       Fastest implementation, using 'hist' and 'reduced.sample'
#
#
# -------- functions ----------------------------------------
#	Kest()		compute estimate of K by border correction
#
#	pairdist()	compute matrix of distances between 	
#			  each pair of data points
#
#       Kount()         internal routine
#
# -------- standard arguments ------------------------------	
#	pp		point pattern (of class 'ppp')
#
#	r		distance values at which to compute K	
#
# -------- standard outputs ------------------------------	
#	border:		K function estimated by border method
#			using standard formula (denominator = count of points)
#
#       bord.modif:	K function estimated by border method
#			using modified formula 
#			(denominator = area of eroded window
#
#	r:		same as input
#
# -------- example -----------------------------------------
#
#	x <- runif(50)
#	y <- runif(50)
#	X <- ppp(x,y,c(0,1),c(0,1))
# either:
#	r <- seq(0,0.25,length=50)
#	k <- Kest(X, r)$border
#	plot(r, k, type="l")
#	lines(r, pi * r^2, lty=2)
# or:
#       Kdat <- Kest(X)
#       plot(border ~ r, type="l", data=Kdat)
# ------------------------------------------------------------------------

"pairdist"<-
function(x, y)
{
	if(length(x) != length(y))
		stop("lengths of x and y do not match")
	n <- length(x)
	xx <- matrix(rep(x, n), nrow = n)
	yy <- matrix(rep(y, n), nrow = n)
	d <- sqrt((xx - t(xx))^2 + (yy - t(yy))^2)
	return(d)
}
	
"Kest"<-
function(X, r=NULL, breaks=NULL, slow=F, both=F, ...)
{
	verifyclass(X, "ppp")

	npoints <- X$n
	area <- area.owin(X$window)
	lambda <- npoints/area
	lambda2 <- (npoints * (npoints - 1))/(area^2)

        breaks <- handle.r.b.args(r, breaks, X$window)
        r <- breaks$r
        
        # pairwise distance
	d <- pairdist(X$x, X$y)
        diag(d) <- Inf

        # distances to boundary
	b <- bdist.points(X)

        # count 
        RS <- Kount(d, b, breaks, slow)
        
        K.b1 <- RS$numerator/(lambda * RS$denom.count)
        
        if(!both)
          return(data.frame(border=K.b1, r=r, theo= pi * r^2))
        
        denom.area <- eroded.areas(X$window, r)
        K.b2 <- RS$numerator/(lambda2 * denom.area)

        return(data.frame(border=K.b1, bord.modif=K.b2, r=r,
                          theo=pi * r^2))
}
	
Kount <- function(d, b, breaks, slow=F) {
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
         nclose <- apply(close, 1, sum) # assumes diag(d) set to Inf
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
        nco <- hist(d[uncen], breaks=breaks$val,plot=F,probability=F)$counts
        # histogram of censoring times for noncensored distances
        ncc <- hist(bb[uncen], breaks=breaks$val,plot=F,probability=F)$counts
        # histogram of censoring times (yes, this is a different total size)
        cen <- hist(b, breaks=breaks$val,plot=F,probability=F)$counts
        # go
        RS <- reduced.sample(nco, cen, ncc, show=T)
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
