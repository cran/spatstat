#
#
#    pairsat.family.S
#
#    $Revision: 1.13 $	$Date: 2003/03/11 07:47:32 $
#
#    The saturated pairwise interaction family of point process models
#
#    (an extension of Geyer's saturation process to all pairwise interactions)
#
#    pairsat.family:         object of class 'isf'
#                     defining saturated pairwise interaction
#	
#
# -------------------------------------------------------------------
#	

pairsat.family <-
  list(
         name  = "saturated pairwise",
         print = function(self) {
                      cat("Saturated pairwise interaction family\n")
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
  #   Equal       logical matrix X[i] == U[j]             matrix or NULL
  #   pairpot     potential function (see above)          function(d, p)
  #   potpars     auxiliary parameters for pairpot        list(......)
  #   correction  edge correction type                    (string)
  #
  #   Note the Geyer saturation threshold must be given in 'potpars$saturate'
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

# coercion should be unnecessary, but this is useful for debugging
X <- as.ppp(X)
U <- as.ppp(U, X$window)   # i.e. X$window is DEFAULT window

# extract coordinates
x <- X$x
y <- X$y
xx <- U$x
yy <- U$y

ndata <- length(x)
nall <- length(xx)
ndummy <- nall - ndata

# saturation parameter
saturate <- potpars$saturate

# validate 'correction' argument
if(!any(correction == c(
          "periodic",
          "border",
          "translate",
          "isotropic", "Ripley",
          "none")))
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

# Sum the pairwise potentials (sum over the data points)
   V <- apply(POT, c(2,3), sum)

if(is.null(saturate))
  return(V)

#################################################################
################## saturation part ##############################
#################################################################

#
# (a) compute SATURATED potential sums
V.sat <- array(pmin(saturate, V), dim=dim(V))

#
# (b) compute effect of addition/deletion of dummy/data point j
# on the UNSATURATED potential sum of each data point i
#
# Identify data points
     # is.data <- apply(Equal, 2, any)
is.data <- matcolany(Equal)              # logical vector corresp. to V

# Extract potential sums for data points only
V.data <- V[is.data, , drop=FALSE]

# replicate them so that V.dat.rep[i,j,k] = V.data[i, k]
V.dat.rep <- aperm(array(V.data, dim=c(dim(V.data), nall)), c(1,3,2))

# make a logical array   col.is.data[i,j,k] = is.data[j]
dip <- dim(POT)
mat <- matrix(, nrow=dip[1], ncol=dip[2])
izdat <- is.data[col(mat)] 
col.is.data <- array(izdat, dim=dip) # automatically replicates

# compute value of unsaturated potential sum for each data point i
# obtained after addition/deletion of each dummy/data point j
                                  
V.after <- V.dat.rep + ifelse(col.is.data, -POT, POT)
#
#
# (c) difference of SATURATED potential sums for each data point i
# before & after increment/decrement of each dummy/data point j
#
# saturated values after increment/decrement
V.after.sat <- array(pmin(saturate, V.after), dim=dim(V.after))
# saturated values before
V.dat.rep.sat <- array(pmin(saturate, V.dat.rep), dim=dim(V.dat.rep))
# difference
V.delta <- V.after.sat - V.dat.rep.sat
V.delta <- ifelse(col.is.data, -V.delta, V.delta)
#
# (d) Sum (c) over all data points i
V.delta.sum <- apply(V.delta, c(2,3), sum)
#
# (e) Result
V <- V.sat + V.delta.sum
  
return(V)

}     ######### end of function $eval                            

)     ######### end of list

class(pairsat.family) <- "isf"
