#
#
#    strauss.S
#
#    $Revision: 2.18 $	$Date: 2010/07/18 08:46:28 $
#
#    The Strauss process
#
#    Strauss()    create an instance of the Strauss process
#                 [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

Strauss <- function(r) {
  out <- 
  list(
       name     = "Strauss process",
       creator  = "Strauss",
       family    = pairwise.family,
       pot      = function(d, par) {
         d <= par$r
       },
       par      = list(r = r),
       parnames = "interaction distance",
       init     = function(self) {
         r <- self$par$r
         if(!is.numeric(r) || length(r) != 1 || r <= 0)
           stop("interaction distance r must be a positive number")
       },
       update = NULL,  # default OK
       print = NULL,    # default OK
       interpret =  function(coeffs, self) {
         loggamma <- as.numeric(coeffs[1])
         gamma <- exp(loggamma)
         return(list(param=list(gamma=gamma),
                     inames="interaction parameter gamma",
                     printable=round(gamma,4)))
       },
       valid = function(coeffs, self) {
         gamma <- ((self$interpret)(coeffs, self))$param$gamma
         return(is.finite(gamma) && (gamma <= 1))
       },
       project = function(coeffs, self) {
         loggamma <- as.numeric(coeffs[1])
         coeffs[1] <- if(is.na(loggamma)) 0 else min(0, loggamma)
         return(coeffs)
       },
       irange = function(self, coeffs=NA, epsilon=0, ...) {
         r <- self$par$r
         if(any(is.na(coeffs)))
           return(r)
         loggamma <- coeffs[1]
         if(abs(loggamma) <= epsilon)
           return(0)
         else
           return(r)
       },
       version=versionstring.spatstat(),
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction, ...) {
         # fast evaluator for Strauss interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Strauss")
         r <- potpars$r
         answer <- strausscounts(U, X, r, EqualPairs)
         return(matrix(answer, ncol=1))
       }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}


strausscounts <- function(U, X, r, EqualPairs=NULL) {
  answer <- crosspaircounts(U,X,r)
  nU <- npoints(U)
  # subtract counts of identical pairs
  if(!is.null(EqualPairs)) {
    idcount <- as.integer(table(factor(EqualPairs[,2], levels=1:nU)))
    answer <- answer - idcount
  }
  return(answer)
}

crosspaircounts <- function(X, Y, r) {
  stopifnot(is.numeric(r))
  # sort in increasing order of x coordinate
  oX <- order(X$x)
  oY <- order(Y$x)
  Xsort <- X[oX]
  Ysort <- Y[oY]
  nX <- npoints(X)
  nY <- npoints(Y)
  # call C routine
  out <- .C("closepaircounts",
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            rrmax    = as.double(r),
            counts   = as.integer(integer(nX)),
            PACKAGE  = "spatstat")
  answer <- integer(nX)
  answer[oX] <- out$counts
  return(answer)
}
