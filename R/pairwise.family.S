#
#
#    pairwise.family.S
#
#    $Revision: 1.8 $	$Date: 2005/04/16 01:35:18 $
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

if(length(correction) > 1)
  stop("Only one edge correction allowed at a time!")

if(!any(correction == c("periodic", "border", "translate", "isotropic", "Ripley", "none")))
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
if(!is.null(Equal) && any(Equal))
  POT[Equal] <- 0

# Sum the pairwise potentials 

V <- apply(POT, c(2,3), sum)

# attach the original pair potentials
attr(V, "POT") <- POT

return(V)

},
######### end of function $eval
       suffstat = function(model, X) {
# for pairwise models only  (possibly nonstationary)
  verifyclass(model, "ppm")
  verifyclass(X, "ppp")

  if(!identical(model$interaction$family$name,"pairwise"))
    stop("Model is not a pairwise interaction process")

  trend <- model$trend
  inter <- model$interaction
  covar <- model$covariates
 
  Empty <- X[rep(FALSE, X$n)]
  Q     <- quad(X, Empty)
  prep  <- mpl.prepare(Q, X, X, trend, inter, covar,
                       correction=model$correction,
                       rbord=model$rbord)
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
  }
######### end of function $suffstat
)
######### end of list

class(pairwise.family) <- "isf"
