#
#
#    multipair.family.S
#
#    $Revision: 1.4 $	$Date: 2002/05/06 09:32:17 $
#
#    Pairwise interaction class for multitype point processes
#    i.e. marked point processes with a finite number of possible types
#
#    multipair.family:      object of class 'isf' defining
#                           pairwise interaction for multitype point processes
#	
# -------------------------------------------------------------------
#	

multipair.family <-
  list(
         name  = "multipair",
         print = function(self) {
                      cat("Multitype pairwise interaction family\n")
         },
         eval  = function(X,U,Equal,pairpot,potpars,correction) {
  #
  # multipair.family$eval
  #
  #  $Revision: 1.4 $  $Date: 2002/05/06 09:32:17 $
  #         
  # This auxiliary function is not meant to be called by the user.
  # It computes the distances between points,
  # evaluates the pair potential and applies edge corrections.
  #
  # Arguments:
  #   X           data point pattern (marked)             'ppp' object
  #   U           points at which to evaluate potential   list(x,y,marks)
  #   Equal       logical matrix X[i] == U[j]             matrix or NULL
  #                      (NB: equality only if marks are equal too)
  #                      NULL means all comparisons are FALSE
  #   pairpot     potential function (see above)          function()
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
  #    pairpot(M, V1, V2, potpars)
  # where M is a matrix of interpoint distances,
  #       V1 and V2 are vectors of marks for the rows and columns of M resp.
  # It must return a matrix with the same dimensions as M
  # or an array with its first two dimensions the same as the dimensions of M.
  ##########################################################################

# coercion should be unnecessary..
# X <- as.ppp(X)
# U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window

x <- X$x
y <- X$y
m <- X$marks

xx <- U$x
yy <- U$y
mm <- U$marks

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

POT <- pairpot(M, m, mm, potpars)
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
        edgewt <- edge.Trans(X, U)
        # sanity check ("everybody knows there ain't no...")
        if(!is.matrix(edgewt))
          stop("internal error: edge.Trans() did not yield a matrix")
        if(nrow(edgewt) != X$n || ncol(edgewt) != length(U$x))
          stop("internal error: edge weights matrix returned by edge.Trans() has wrong dimensions")
        POT <- c(edgewt) * POT
} else if(correction == "isotropic" || correction == "Ripley") {
        # weights are required for contributions from QUADRATURE points
        edgewt <- t(edge.Ripley(U, t(M), X$window))
        if(!is.matrix(edgewt))
          stop("internal error: edge.Ripley() did not return a matrix")
        if(nrow(edgewt) != X$n || ncol(edgewt) != length(U$x))
          stop("internal error: edge weights matrix returned by edge.Ripley() has wrong dimensions")
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

class(multipair.family) <- "isf"





###########    utilities for this family  #######################


MultiPair.checkmatrix <-
  function(mat, n, name) {
    if(!is.matrix(mat))
      stop(paste(name, "must be a matrix"))
    if(any(dim(mat) != rep(n,2)))
      stop(paste(name, "must be a square matrix,",
                 "of size", n, "x", n))
    isna <- is.na(mat)
    if(any(mat[!isna] <= 0))
      stop(paste("Entries of", name,
                 "must be positive numbers or NA"))
    if(any(isna != t(isna)) ||
       any(mat[!isna] != t(mat)[!isna]))
      stop(paste(name, "must be a symmetric matrix"))
  }

