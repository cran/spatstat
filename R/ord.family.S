#
#
#    ord.family.S
#
#    $Revision: 1.8 $	$Date: 2001/08/09 10:53:16 $
#
#    The Ord model (family of point process models)
#
#    ord.family:      object of class 'isf' defining Ord model structure
#	
#
# -------------------------------------------------------------------
#	

ord.family <-
  list(
         name  = "ord",
         print = function(self) {
                      cat("Ord model family\n")
         },
         eval  = function(X, U, Equal, pot, pars, ...) {
  #
  # This auxiliary function is not meant to be called by the user.
  # It computes the distances between points,
  # evaluates the pair potential and applies edge corrections.
  #
  # Arguments:
  #   X           data point pattern                      'ppp' object
  #   U           points at which to evaluate potential   list(x,y) suffices
  #   Equal       logical matrix X[i] == U[j]             matrix or NULL
  #   pot         potential function                      function(d, p)
  #   pars        auxiliary parameters for pot            list(......)
  #   ...         IGNORED                             
  #
  # Value:
  #    matrix of values of the potential
  #    induced by the pattern X at each location given in U.
  #    The rows of this matrix correspond to the rows of U (the sample points);
  #    the k columns are the coordinates of the k-dimensional potential.
  #
  # Note:
  # The potential function 'pot' will be called as
  #    pot(M, pars)   where M is a vector of tile areas.
  # It must return a vector of the same length as M
  # or a matrix with number of rows equal to the length of M
  ##########################################################################

nall <- length(U$x)       # number of data + dummy points

# determine which points in the combined list are data points
if(!is.null(Equal))           
  #is.data <- apply(Equal, 2, any)
  is.data <- matcolany(Equal)
else
  is.data <- rep(F, nall)

#############################################################################
# First compute Dirichlet tessellation of data
# and its total potential (which could be vector-valued)
#############################################################################

Wdata <- dirichlet.weights(X)   # sic - these are the tile areas.
Pdata <- pot(Wdata, pars)
summa <- function(P) {
  if(is.vector(P)) sum(P) else matrowsum(P)
}
total.data.potential <- summa(Pdata)

# Initialise V

dimpot <- dim(Pdata)[-1]  # dimension of each value of the potential function
                          # (NULL if scalar)
dimV <- c(nall, dimpot)
V <- array(0, dim=dimV)

rowV <- array(1:nall, dim=dimV)

#################### Next, evaluate V for the data points.  ###############
# For each data point, compute Dirichlet tessellation
# of the data with this point removed.
# Compute difference of total potential.
#############################################################################


Xminus <- X
for(j in seq(X$n)) {
  Xminus$x <- X$x[-j]
  Xminus$y <- X$y[-j]
        #  Dirichlet tessellation of data without point j
  Wminus <- dirichlet.weights(Xminus)
        #  regressor is the difference in total potential
  V[rowV == j] <- total.data.potential - summa(pot(Wminus, pars))
}


#################### Next, evaluate V for the dummy points   ################
# For each dummy point, compute Dirichlet tessellation
# of (data points together with this dummy point) only. 
# Take difference of total potential.
#############################################################################

Xplus <- X

for(j in seq(U$x)[!is.data]) {
		Xplus$x <- c(U$x[j], X$x)
		Xplus$y <- c(U$y[j], X$y)
                     #  compute Dirichlet tessellation (of these points only!)
                Wplus <- dirichlet.weights(Xplus)
                     #  regressor is difference in total potential
                V[rowV == j] <- summa(pot(Wplus)) - total.data.potential
}

return(V)

} ######### end of function $eval                            

) ######### end of list

class(ord.family) <- "isf"
