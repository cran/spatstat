#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 5.2 $	$Date: 2002/05/13 12:41:10 $
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
#	crossdist()	compute matrix of distances between 	
#			  each pair of data points 
#			  in two separate lists of points
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

"crossdist"<-
function(x1, y1, x2, y2)
{
        # returns matrix[i,j] = distance from (x1[i],y1[i]) to (x2[j],y2[j])
	if(length(x1) != length(y1))
		stop("lengths of x1 and y1 do not match")
	if(length(x2) != length(y2))
		stop("lengths of x2 and y2 do not match")
	n1 <- length(x1)
	n2 <- length(x2)
	X1 <- matrix(rep(x1, n2), ncol = n2)
	Y1 <- matrix(rep(y1, n2), ncol = n2)
	X2 <- matrix(rep(x2, n1), ncol = n1)
	Y2 <- matrix(rep(y2, n1), ncol = n1)
	d <- sqrt((X1 - t(X2))^2 + (Y1 - t(Y2))^2)
	return(d)
}

"Kcross" <- 
function(X, i=1, j=2, r=NULL, breaks=NULL,
         correction =c("border", "isotropic", "Ripley", "translate") , ...)
{
    
	verifyclass(X, "ppp")
	if(!is.marked(X))
		stop("point pattern has no marks (no component 'marks')")
	
	I <- (X$marks == i)
	J <- (X$marks == j)
	
	if(!any(I)) stop(paste("No points have mark i =", i))
	if(!any(J)) stop(paste("No points have mark j =", j))
	
	Kmulti(X, I, J, r, breaks)
}

"Kdot" <- 
function(X, i=1, r=NULL, breaks=NULL,
         correction = c("border", "isotropic", "Ripley", "translate") , ...)
{
	verifyclass(X, "ppp")
	if(!is.marked(X))
		stop("point pattern has no marks (no component 'marks')")
	
	I <- (X$marks == i)
	J <- rep(TRUE, X$n)  # i.e. all points
	
	if(!any(I)) stop(paste("No points have mark i =", i))
	
	Kmulti(X, I, J, r, breaks)
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

        breaks <- handle.r.b.args(r, breaks, W)
        r <- breaks$r
        
        # available selection of edge corrections depends on window
        if(W$type != "rectangle") {
           iso <- (correction == "isotropic") | (correction == "Ripley")
           if(all(iso))
             stop("Isotropic correction not implemented for non-rectangular windows")
           if(any(iso)) {
             if(!missing(correction))
               warning("Isotropic correction not implemented for non-rectangular windows")
             correction <- correction[!iso]
           }
        }
         
	if(!is.logical(I) || !is.logical(J))
		stop("I and J must be logical vectors")
	if(length(I) != npoints || length(J) != npoints)
	     stop("The length of I and J must equal \
 the number of points in the pattern")
	
	if(!any(I)) stop("no points satisfy I")
	if(!any(J)) stop("no points satisfy J")
		
	nI <- sum(I)
	nJ <- sum(J)
	lambdaI <- nI/area
	lambdaJ <- nJ/area

# This will be the output data frame
# It will be given more columns later        
        K <- data.frame(r=r, theo= pi * r^2)
        
# interpoint distances		
	d <- crossdist(x[I], y[I], x[J], y[J])
# distances to boundary	
	b <- (bdist.points(X))[I]
        
# Determine which interpoint distances d[i,j] refer to the same point
# (not just which distances are zero)        
        identical <- matrix(FALSE, nrow=nI, ncol=nJ)
        common <- I & J
        if(any(common)) {
          Irow <- cumsum(I)
          Jcol <- cumsum(J)
          icommon <- (1:npoints)[common]
          for(i in icommon)
            identical[Irow[i], Jcol[i]] <- TRUE
        }

# Compute estimates by each of the selected edge corrections.
        
        if(any(correction == "border" | correction == "bord.modif")) {
          # border method
          # Compute distances to boundary
          b <- bdist.points(X[I])
          # Distances corresponding to identical pairs
          # are excluded from consideration
          d[identical] <- Inf
          # apply reduced sample algorithm
          RS <- Kount(d, b, breaks, slow=FALSE)
          if(any(correction == "border"))
            K$border <- RS$numerator/(lambdaJ * RS$denom.count)
          if(any(correction == "bord.modif")) {
            denom.area <- eroded.areas(W, r)
            K$bord.modif <- RS$numerator/(denom.area * nI * nJ)
          }
          # reset identical pairs to original values
          d[identical] <- 0
        }
        if(any(correction == "translate")) {
          # translation correction
            edgewt <- edge.Trans(X[I], X[J])
            wh <- whist(d[!identical], breaks$val, edgewt[!identical])
            K$trans <- cumsum(wh)/(lambdaI * lambdaJ * area)
            rmax <- diameter(W)/2
            K$trans[r >= rmax] <- NA
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(X[I], d)
            wh <- whist(d[!identical], breaks$val, edgewt[!identical])
            K$iso <- cumsum(wh)/(lambdaI * lambdaJ * area)
            rmax <- diameter(W)/2
            K$iso[r >= rmax] <- NA
        }

        return(K)
          
}
