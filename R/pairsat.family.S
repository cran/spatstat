#
#
#    pairsat.family.S
#
#    $Revision: 1.26 $	$Date: 2007/03/16 08:35:55 $
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
         eval  = function(X,U,Equal,pairpot,potpars,correction,
                          ..., halfway=FALSE) {
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
  #   Note the Geyer saturation threshold must be given in 'potpars$sat'
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

# saturation parameter(s)
saturate <- potpars$sat

if(is.null(saturate)) {
  # pairwise interaction 
  V <- pairwise.family$eval(X, U, Equal, pairpot, potpars, correction, ...)
  return(V)
}

# first check that all data points are included in the quadrature points
missingdata <- !as.vector(matrowany(Equal))
somemissing <- any(missingdata)
if(somemissing) {
  # add the missing data points
  U <- superimpose(U, X[missingdata], W=X$window)
  # extend the equality matrix to maintain Equal[i,j] = (X[i] == U[j])
  originalcolumns <- seq(ncol(Equal))
  sn <- seq(X$n)
  id <- outer(sn, sn[missingdata], "==")
  Equal <- cbind(Equal, id)
}

# compute the pair potentials POT and the unsaturated potential sums V

V <- pairwise.family$eval(X, U, Equal, pairpot, potpars, correction, ...)
POT <- attr(V, "POT")

computed <- attr(V, "computed")   # could be NULL

#
# V is a matrix with rows = quadrature points,
#                    columns = coordinates of potential
# POT is an array with rows = data points
#                      columns = quadrature points
#                      planes = coordinates of potential

#################################################################
################## saturation part ##############################
#################################################################

# check dimensions and ensure 'saturate' is a vector
ns <- length(saturate)
np <- ncol(V)
if(ns == 1 && np > 1)
  saturate <- rep(saturate, np)
else if(ns != np)
  stop("Length of vector of saturation parameters is incompatible with the pair potential", call.=FALSE)

# replicate as a matrix and as an array
saturate2 <- array(saturate[slice.index(V, 2)], dim=dim(V))
saturate3 <- array(saturate[slice.index(POT, 3)], dim=dim(POT))
#
# (a) compute SATURATED potential sums
V.sat <- pmin(V, saturate2)

if(halfway)
  return(V.sat)
#
# (b) compute effect of addition/deletion of dummy/data point j
# on the UNSATURATED potential sum of each data point i
#
# Identify data points
     # is.data <- apply(Equal, 2, any)
is.data <- as.vector(matcolany(Equal))  # logical vector corresp. to rows of V

# Extract potential sums for data points only
V.data <- V[is.data, , drop=FALSE]

# replicate them so that V.dat.rep[i,j,k] = V.data[i, k]
V.dat.rep <- aperm(array(V.data, dim=c(dim(V.data), U$n)), c(1,3,2))

# make a logical array   col.is.data[i,j,k] = is.data[j]
col.is.data <- array(is.data[slice.index(POT, 2)], dim=dim(POT))

# compute value of unsaturated potential sum for each data point i
# obtained after addition/deletion of each dummy/data point j
                                  
V.after <- V.dat.rep + ifelse(col.is.data, -POT, POT)
#
#
# (c) difference of SATURATED potential sums for each data point i
# before & after increment/decrement of each dummy/data point j
#
# saturated values after increment/decrement
V.after.sat <- array(pmin(saturate3, V.after), dim=dim(V.after))
# saturated values before
V.dat.rep.sat <- array(pmin(saturate3, V.dat.rep), dim=dim(V.dat.rep))
# difference
V.delta <- V.after.sat - V.dat.rep.sat
V.delta <- ifelse(col.is.data, -V.delta, V.delta)
#
# (d) Sum (c) over all data points i
V.delta.sum <- apply(V.delta, c(2,3), sum)
#
# (e) Result
V <- V.sat + V.delta.sum

##########################################
# remove any columns that were added
if(somemissing)
      V <- V[originalcolumns, , drop=FALSE]

### tack on the saved computations from pairwise.family$eval
attr(V, "computed") <- computed

return(V)

},     ######### end of function $eval                            
suffstat = function(model, X, callstring="pairsat.family$suffstat") {

# for saturated pairwise models only  (possibly nonstationary)
  verifyclass(model, "ppm")
  verifyclass(X, "ppp")

  if(!identical(model$interaction$family$name,"saturated pairwise"))
    stop("Model is not a saturated pairwise interaction process")

  trend <- model$trend
  inter <- model$interaction
  covar <- model$covariates
 
  Empty <- X[rep(FALSE, X$n)]
  Q     <- quad(X, Empty)
  prep  <- mpl.prepare(Q, X, X, trend, inter, covar,
                       correction=model$correction,
                       rbord=model$rbord, Pname="data points",
                       callstring=callstring,
                       halfway=TRUE)
  # halfway=TRUE is passed to inter$family$eval
  # and yields matrix of saturated potential sums 
  fmla    <- prep$fmla
  glmdata <- prep$glmdata

  mof <- model.frame(fmla, glmdata)
  mom <- model.matrix(fmla, mof)

  if(any(colnames(mom) != names(coef(model))))
    warning("Internal error: mismatch between column names of model matrix and names of coefficient vector in fitted model")
     
  order2  <- colnames(mom) %in% model$internal$Vnames
  order1  <- !order2

  if(any(order1)) {
    o1terms <- mom[, order1, drop=FALSE]
    o1sum   <- apply(o1terms, 2, sum)
  }
  if(any(order2)) {
    o2terms <- mom[, order2, drop=FALSE]
    o2sum   <- apply(o2terms, 2, sum)/2
  }

  result <- 0 * coef(model)
  result[order1] <- o1sum
  result[order2] <- o2sum
  return(result)
         

} ######### end of function $suffstat
)     ######### end of list

class(pairsat.family) <- "isf"
