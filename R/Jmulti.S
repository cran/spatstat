#	Jmulti.S
#
#	Usual invocations to compute multitype J function(s)
#	if F and G are not required 
#
#	$Revision: 4.3 $	$Date: 2001/11/12 10:40:04 $
#
#
#
"Jcross" <-
function(X, i=1, j=2, eps=NULL, r=NULL, breaks=NULL) {
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
          stop("point pattern has no \'marks\'")
        I <- (X$marks == i)
        J <- (X$marks == j)
        result <- Jmulti(X, I, J, eps, r, breaks)
	return(result)
}

"Jdot" <-
function(X, i=1, eps=NULL, r=NULL, breaks=NULL) {
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
          stop("point pattern has no \'marks\'")
        I <- (X$marks == i)
        J <- rep(T, X$n)
        result <- Jmulti(X, I, J, eps, r, breaks)
	return(result)
}

"Jmulti" <- 	
function(X, I, J, eps=NULL, r=NULL, breaks=NULL) {
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
        brks <- handle.r.b.args(r, breaks, X$window, eps)$val
	FJ <- Fest(X[J], eps, breaks=brks)
	GIJ <- Gmulti(X, I, J, breaks=brks)
        ratio <- function(a, b, c) {
          result <- a/b
          result[ b == 0 ] <- c
          result
        }
	Jkm <- ratio(1-GIJ$km, 1-FJ$km, NA)
	Jrs <- ratio(1-GIJ$rs, 1-FJ$rs, NA)
	Jun <- ratio(1-GIJ$raw, 1-FJ$raw, NA)
        theo <- rep(1, length(FJ$r))
        
	result <- data.frame(J=Jkm,r=FJ$r,km=Jkm,rs=Jrs,un=Jun,theo=theo)
        attr(result, "G") <- GIJ
        attr(result, "F") <- FJ
        result
}
