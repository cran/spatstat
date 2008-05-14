#	Jmulti.S
#
#	Usual invocations to compute multitype J function(s)
#	if F and G are not required 
#
#	$Revision: 4.18 $	$Date: 2008/04/28 23:21:05 $
#
#
#
"Jcross" <-
function(X, i=1, j=2, eps=NULL, r=NULL, breaks=NULL, ...) {
#
#       multitype J function J_{ij}(r)
#  
#	X:		point pattern (an object of class 'ppp')
#       i, j:           types for which J_{i,j}(r) is calculated  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
        X <- as.ppp(X)
        if(!is.marked(X))
          stop(paste("point pattern has no", sQuote("marks")))

        marx <- marks(X)
        I <- (marx == i)
        if(sum(I) == 0)
          stop(paste("No points have mark = ", i))
        
        if(i == j)
          result <- Jest(X[I], eps=eps, r=r, breaks=breaks)
        else {
          J <- (marx == j)
          result <- Jmulti(X, I, J,
                           eps=eps, r=r, breaks=breaks, disjoint=TRUE)
        }
        attr(result, "ylab") <- substitute(Jcross[i,j](r), list(i=paste(i),j=paste(j)))
	return(result)
}

"Jdot" <-
function(X, i=1, eps=NULL, r=NULL, breaks=NULL, ...) {
#  
#    multitype J function J_{i\dot}(r)
#  
#	X:		point pattern (an object of class 'ppp')
#       i:              mark i for which we calculate J_{i\cdot}(r)  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#
        X <- as.ppp(X)
        if(!is.marked(X))
          stop(paste("point pattern has no", sQuote("marks")))
        I <- (marks(X) == i)
        if(sum(I) == 0)
          stop(paste("No points have mark = ", i))          
        J <- rep(TRUE, X$n)
        result <- Jmulti(X, I, J,
                         eps=eps, r=r, breaks=breaks, disjoint=FALSE)
        attr(result, "ylab") <- substitute(Jdot[i](r), list(i=paste(i)))
	return(result)
}

"Jmulti" <- 	
function(X, I, J, eps=NULL, r=NULL, breaks=NULL, ..., disjoint=NULL) {
#  
#    multitype J function (generic engine)
#  
#	X		marked point pattern (of class ppp)
#	
#	I,J		logical vectors of length equal to the number of points
#			and identifying the two subsets of points to be
#			compared.
#  
#	eps:		raster grid mesh size for distance transform
#				(unless specified by X$window)
#  
#       r:              (optional) values of argument r  
#	breaks:		(optional) breakpoints for argument r
#  
#
        X <- as.ppp(X)
        W<- X$window
        rmaxdefault <- rmax.rule("J", W)
        brks <- handle.r.b.args(r, breaks, W, rmaxdefault=rmaxdefault)$val
	FJ <- Fest(X[J], eps, breaks=brks)
	GIJ <- Gmulti(X, I, J, breaks=brks, disjoint=disjoint)
        ratio <- function(a, b, c) {
          result <- a/b
          result[ b == 0 ] <- c
          result
        }
	Jkm <- ratio(1-GIJ$km, 1-FJ$km, NA)
	Jrs <- ratio(1-GIJ$rs, 1-FJ$rs, NA)
	Jun <- ratio(1-GIJ$raw, 1-FJ$raw, NA)
        theo <- rep(1, length(FJ$r))
        
	result <- data.frame(r=FJ$r,theo=theo,un=Jun,rs=Jrs,km=Jkm)

        alim <- range(result$r[FJ$km <= 0.9])
        labl <- c("r", "Jpois(r)", "Jun(r)", "Jbord(r)", "Jkm(r)")
        desc <- c("distance argument r",
                  "theoretical Poisson J(r)=1",
                  "uncorrected estimate of J(r)",
                  "border corrected estimate of J(r)",
                  "Kaplan-Meier estimate of J(r)")
        Z <- fv(result,
                "r", substitute(Jmulti(r), NULL), "km", . ~ r, alim, labl, desc)
        
        attr(Z, "G") <- GIJ
        attr(Z, "F") <- FJ
        unitname(Z) <- unitname(X)
        return(Z)
}
