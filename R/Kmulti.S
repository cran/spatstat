#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 4.4 $	$Date: 2001/11/06 11:05:09 $
#
#	Border correction only.
#	Arbitrary window.
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
# -------- standard outputs ------------------------------	
#	border:		cross K function estimated by border method
#			using standard formula (denominator = count of points)
#
#       bord.modif:	cross K function estimated by border method
#			using modified formula 
#			(denominator = area of eroded window
#
#	r:		same as input
#
#
# -------- example -----------------------------------------
#
#	x <- runif(50)
#	y <- runif(50)
#	m <- c(rep(1,25),rep(2,25))
#	X <- ppp(x,y,c(0,1),c(0,1),marks=m)
# either:
#	r <- seq(0,0.5,length=50)
#	k <- Kcross(X, 1, 2, r)$border
#	plot(r, k, type="l")
#	lines(r, pi * r^2, lty=2)
# or:
#       Kdat <- Kcross(X, 1, 2)
#       plot(border ~ r, type="l", data=Kdat)
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
function(X, i=1, j=2, r=NULL, breaks=NULL, ...) {
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
function(X, i=1, r=NULL, breaks=NULL, ...) {
	verifyclass(X, "ppp")
	if(!is.marked(X))
		stop("point pattern has no marks (no component 'marks')")
	
	I <- (X$marks == i)
	J <- rep(T, X$n)  # i.e. all points
	
	if(!any(I)) stop(paste("No points have mark i =", i))
	
	Kmulti(X, I, J, r, breaks)
}


"Kmulti"<-
function(X, I, J, r=NULL, breaks=NULL, ...)
{
	verifyclass(X, "ppp")

	npoints <- X$n
	x <- X$x
	y <- X$y
	area <- area.owin(X$window)

        breaks <- handle.r.b.args(r, breaks, X$window)
        r <- breaks$r
        
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

# interpoint distances		
	d <- crossdist(x[I], y[I], x[J], y[J])
# distances to boundary	
	b <- (bdist.points(X))[I]
        
# interpoint distances that refer to identical points
# are excluded from consideration by setting them to infinity        
        common <- I & J
        if(any(common)) {
          Irow <- cumsum(I)
          Jcol <- cumsum(J)
          icommon <- (1:npoints)[common]
          for(i in icommon)
            d[Irow[i], Jcol[i]] <- Inf
        }
          
# border correction estimator          
	RS <- Kount(d, b, breaks, slow=F)  
        
        K.b1 <- RS$numerator/(lambdaJ * RS$denom.count)

# variant 
        denom.area <- eroded.areas(X$window, r)
        K.b2 <- RS$numerator/(denom.area * nI * nJ)

# theoretical value for marked Poisson
        theo <- pi * r^2
        
        return(data.frame(border=K.b1, bord.modif = K.b2, r=r, theo=theo))
          
}
