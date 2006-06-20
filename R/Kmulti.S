#
#	Kmulti.S		
#
#	Compute estimates of cross-type K functions
#	for multitype point patterns
#
#	$Revision: 5.14 $	$Date: 2006/06/11 23:10:21 $
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
	
	result <- Kmulti(X, I, J, r, breaks)
        attr(result, "ylab") <- substitute(Kcross(r), NULL)
        result
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
	
	result <- Kmulti(X, I, J, r, breaks)
        attr(result, "ylab") <- substitute(Kdot(r), NULL)
        result
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

        # available selection of edge corrections depends on window
        if(W$type == "mask") {
           iso <- (correction == "isotropic") | (correction == "Ripley")
           if(all(iso))
             stop("Isotropic correction not implemented for binary masks")
           if(any(iso)) {
             if(!missing(correction))
               warning("Isotropic correction not implemented for binary masks")
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
        desc <- c("distance argument r", "theoretical Poisson K(r)")
        K <- fv(K, "r", substitute(Kmulti(r), NULL),
                "theo", , alim, c("r","Kpois(r)"), desc)

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
            K <- bind.fv(K, data.frame(bord.modif=Kbm), "Kbord*(r)",
                         "modified border-corrected estimate of K(r)",
                         "bord.modif")
          }
          if(any(correction == "border")) {
            Kb <- RS$numerator/(lambdaJ * RS$denom.count)
            K <- bind.fv(K, data.frame(border=Kb), "Kbord(r)",
                         "border-corrected estimate of K(r)",
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
            K <- bind.fv(K, data.frame(trans=Ktrans), "Ktrans(r)",
                         "translation-corrected estimate of K(r)",
                         "trans")
        }
        if(any(correction == "isotropic" | correction == "Ripley")) {
          # Ripley isotropic correction
            edgewt <- edge.Ripley(XI[icloseI], matrix(dcloseIJ, ncol=1))
            wh <- whist(dcloseIJ, breaks$val, edgewt)
            Kiso <- cumsum(wh)/(lambdaI * lambdaJ * area)
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
