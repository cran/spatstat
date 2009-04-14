#	Jmulti.S
#
#	Usual invocations to compute multitype J function(s)
#	if F and G are not required 
#
#	$Revision: 4.22 $	$Date: 2009/04/14 23:31:15 $
#
#
#
"Jcross" <-
function(X, i, j, eps=NULL, r=NULL, breaks=NULL, ...) {
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
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
  if(missing(j)) j <- levels(marx)[2]
#
  I <- (marx == i)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))
#        
  if(i == j)
    result <- Jest(X[I], eps=eps, r=r, breaks=breaks)
  else {
    J <- (marx == j)
    result <- Jmulti(X, I, J,
                     eps=eps, r=r, breaks=breaks, disjoint=TRUE)
  }
  result <-
    rebadge.fv(result,
               substitute(Jcross[i,j](r), list(i=paste(i),j=paste(j))),
               "J*",
               new.yexp=substitute(Jcross[list(i,j)](r),
                                   list(i=paste(i),j=paste(j))))
  return(result)
}

"Jdot" <-
function(X, i, eps=NULL, r=NULL, breaks=NULL, ...) {
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
  stopifnot(is.multitype(X))
#
  marx <- marks(X, dfok=FALSE)
  if(missing(i)) i <- levels(marx)[1]
#  
  I <- (marx == i)
  if(sum(I) == 0)
    stop(paste("No points have mark = ", i))          
  J <- rep(TRUE, X$n)
#  
  result <- Jmulti(X, I, J,
                   eps=eps, r=r, breaks=breaks, disjoint=FALSE)
  result <-
    rebadge.fv(result,
            substitute(Jdot[i](r), list(i=paste(i))),
            "J.")
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
  labl <- c("r", "%spois(r)", "%sun(r)", "%sbord(r)", "%skm(r)")
  desc <- c("distance argument r",
            "theoretical Poisson %s",
            "uncorrected estimate of %s",
            "border corrected estimate of %s",
            "Kaplan-Meier estimate of %s")
  Z <- fv(result,
          "r", substitute(Jmulti(r), NULL), "km", . ~ r, alim, labl, desc,
          fname="Jm")
        
  attr(Z, "G") <- GIJ
  attr(Z, "F") <- FJ
  unitname(Z) <- unitname(X)
  return(Z)
}
