#
#     dg.S
#
#    $Revision: 1.11 $	$Date: 2010/07/18 08:46:28 $
#
#     Diggle-Gratton pair potential
#
#
DiggleGratton <- function(delta, rho) {
  out <- 
  list(
         name     = "Diggle-Gratton process",
         creator  = "DiggleGratton",
         family    = pairwise.family,
         pot      = function(d, par) {
                       delta <- par$delta
                       rho <- par$rho
                       above <- (d > rho)
                       inrange <- (!above) & (d > delta)
                       h <- above + inrange * (d - delta)/(rho - delta)
                       return(log(h))
                    },
         par      = list(delta=delta, rho=rho),
         parnames = list("lower limit delta", "upper limit rho"),
         init     = function(self) {
                      r <- self$par$delta
                      r <- self$par$rho
                      if(!is.numeric(delta) || length(delta) != 1)
                       stop("lower limit delta must be a single number")
                      if(!is.numeric(rho) || length(rho) != 1)
                       stop("upper limit rho must be a single number")
                      stopifnot(delta >= 0)
                      stopifnot(rho > delta)
                      stopifnot(is.finite(rho))
                    },
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           kappa <- as.numeric(coeffs[1])
           return(list(param=list(kappa=kappa),
                       inames="exponent kappa",
                       printable=round(kappa,4)))
         },
         valid = function(coeffs, self) {
           kappa <- ((self$interpret)(coeffs, self))$param$kappa
           return(is.finite(kappa) && (kappa >= 0))
         },
         project = function(coeffs, self) {
           kappa <- as.numeric(coeffs[1])
           coeffs[1] <- if(is.na(kappa)) 0 else max(0, kappa)
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           rho <- self$par$rho
           if(all(is.na(coeffs)))
             return(rho)
           kappa <- coeffs[1]
           delta <- self$par$delta
           if(abs(kappa) <= epsilon)
             return(delta)
           else return(rho)
         },
       version=versionstring.spatstat(),
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction, ...) {
         # fast evaluator for DiggleGratton interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for DiggleGratton")
         delta <- potpars$delta
         rho   <- potpars$rho
         idX <- seq(npoints(X))
         idU <- rep(-1, npoints(U))
         idU[EqualPairs[,2]] <- EqualPairs[,1]
         answer <- diggraterms(U, X, idU, idX, delta, rho)
         answer <- log(pmax(0, answer))
         return(matrix(answer, ncol=1))
       }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}


diggraterms <- function(X, Y, idX, idY, delta, rho) {
  stopifnot(is.numeric(delta))
  stopifnot(is.numeric(rho))
  stopifnot(delta < rho)
  # sort in increasing order of x coordinate
  oX <- order(X$x)
  oY <- order(Y$x)
  Xsort <- X[oX]
  Ysort <- Y[oY]
  idXsort <- idX[oX]
  idYsort <- idY[oY]
  nX <- npoints(X)
  nY <- npoints(Y)
  # call C routine
  out <- .C("Ediggra",
            nnsource = as.integer(nX),
            xsource  = as.double(Xsort$x),
            ysource  = as.double(Xsort$y),
            idsource = as.integer(idXsort),
            nntarget = as.integer(nY),
            xtarget  = as.double(Ysort$x),
            ytarget  = as.double(Ysort$y),
            idtarget = as.integer(idYsort),
            ddelta   = as.double(delta),
            rrho     = as.double(rho),
            values   = as.double(double(nX)),
            PACKAGE  = "spatstat")
  answer <- integer(nX)
  answer[oX] <- out$values
  return(answer)
}
