#
#
#    pairwise.family.S
#
#    $Revision: 1.4 $	$Date: 2001/08/07 08:44:48 $
#
#    The pairwise interaction family of point process models
#
#    pairwise.family:      object of class 'isf' defining pairwise interaction
#	
#
# -------------------------------------------------------------------
#	

pairwise.family <-
  list(
         name  = "pairwise",
         print = function(self) {
                      cat("Pairwise interaction family\n")
         },
         eval  = function(X,U,Equal,pairpot,potpars,correction) {
  #
  # This auxiliary function is not meant to be called by the user.
  # It computes the distances between points,
  # evaluates the pair potential and applies edge corrections.
  #
  # Arguments:
  #   X           data point pattern                      'ppp' object
  #   U           points at which to evaluate potential   list(x,y) suffices
  #   Equal        logical matrix X[i] == U[j]             matrix or NULL
  #   pairpot     potential function (see above)          function(d, p)
  #   potpars     auxiliary parameters for pairpot        list(......)
  #   correction  edge correction type                    (string)
  #
  # Value:
  #    matrix of values of the total pair potential
  #    induced by the pattern X at each location given in U.
  #    The rows of this matrix correspond to the rows of U (the sample points);
  #    the k columns are the coordinates of the k-dimensional potential.
  #
  # Note:
  # The pair potential function 'pairpot' will be called as
  #    pairpot(M, potpars)   where M is a matrix of interpoint distances.
  # It must return a matrix with the same dimensions as M
  # or an array with its first two dimensions the same as the dimensions of M.
  ##########################################################################

# coercion should be unnecessary..
# X <- as.ppp(X)
# U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window

x <- X$x
y <- X$y
xx <- U$x
yy <- U$y
  
if(!any(correction == c("periodic", "border", "translate", "none")))
  stop(paste("Unrecognised edge correction \'", correction, "\'", sep=""))
          
#  
# Form the matrix of distances
	
sqdif <- function(u,v) {(u-v)^2}

MX <- outer(x,xx,sqdif)
MY <- outer(y,yy,sqdif)

if(correction=="periodic") {
	if(X$window$type != "rectangle")
		stop("Periodic edge correction can't be applied",
                     "in an irregular window")
	wide <- diff(X$window$xrange)
	high <- diff(X$window$yrange)
	MX1 <- outer(x,xx-wide,sqdif)
	MX2 <- outer(x,xx+wide,sqdif)
	MX <- pmin(MX, MX1, MX2)
	MY1 <- outer(y,yy-high,sqdif)
	MY2 <- outer(y,yy+high,sqdif)
	MY <- pmin(MY, MY1, MY2)
}
M <- sqrt(MX + MY)

# Evaluate the pairwise potential 

POT <- pairpot(M, potpars)
if(length(dim(POT)) == 1 || any(dim(POT)[1:2] != dim(M))) {
        whinge <- paste(
           "The pair potential function ",deparse(substitute(pairpot)),
           "must produce a matrix or array with its first two dimensions\n",
           "the same as the dimensions of its input.\n", sep="")
	stop(whinge)
}

# make it a 3D array
if(length(dim(POT))==2)
        POT <- array(POT, dim=c(dim(POT),1), dimnames=NULL)
                          
if(correction == "translate") {
	if(X$window$type != "rectangle")
		stop("sorry, translation correction is not yet implemented",
                     "for irregular windows")
	wide <- diff(X$window$xrange)
	high <- diff(X$window$yrange)
        DX <- abs(outer(x,xx,"-"))
        DY <- abs(outer(y,yy,"-"))
                                        # translation correction
        edgewt <- wide * high / ((wide - DX) * (high - DY))
        edgewt <- pmin(edgewt, 1000)    # arbitrary stabilisation
        POT <- c(edgewt) * POT
}

# No pair potential term between a point and itself
if(!is.null(Equal))
  POT[Equal] <- 0

# Sum the pairwise potentials 

V <- apply(POT, c(2,3), sum)

return(V)

}
######### end of function $eval                            
)
######### end of list

class(pairwise.family) <- "isf"
