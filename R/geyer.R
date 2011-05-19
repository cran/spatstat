#
#
#    geyer.S
#
#    $Revision: 2.17 $	$Date: 2011/05/18 02:06:52 $
#
#    Geyer's saturation process
#
#    Geyer()    create an instance of Geyer's saturation process
#                 [an object of class 'interact']
#
#	

Geyer <- function(r, sat) {
  out <- 
  list(
         name     = "Geyer saturation process",
         creator  = "Geyer",
         family    = pairsat.family,
         pot      = function(d, par) {
                         (d <= par$r)  # same as for Strauss
                    },
         par      = list(r = r, sat=sat),
         parnames = c("interaction distance","saturation parameter"),
         init     = function(self) {
                      r <- self$par$r
                      sat <- self$par$sat
                      if(!is.numeric(r) || length(r) != 1 || r <= 0)
                       stop("interaction distance r must be a positive number")
                      if(!is.numeric(sat) || length(sat) != 1 || sat < 1)
                       stop("saturation parameter sat must be a number >= 1")
                      if(ceiling(sat) != floor(sat))
                        warning(paste("saturation parameter sat",
                                      "has a non-integer value"))
                    },
         update = NULL, # default OK
         print = NULL,    # default OK
         interpret =  function(coeffs, self) {
           loggamma <- as.numeric(coeffs[1])
           gamma <- exp(loggamma)
           return(list(param=list(gamma=gamma),
                       inames="interaction parameter gamma",
                       printable=round(gamma,4)))
         },
         valid = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           return(is.finite(gamma))
         },
         project = function(coeffs, self) {
           gamma <- (self$interpret)(coeffs, self)$param$gamma
           if(is.na(gamma)) 
             coeffs[1] <- 0
           else if(!is.finite(gamma)) 
             coeffs[1] <- log(.Machine$double.xmax)/self$par$sat
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$r
           if(any(!is.na(coeffs))) {
             loggamma <- coeffs[1]
             if(!is.na(loggamma) && (abs(loggamma) <= epsilon))
               return(0)
           }
           return(2 * r)
         },
       version=versionstring.spatstat(),
       # fast evaluation is available for the border correction only
       can.do.fast=function(X,correction,par) {
         return(all(correction %in% c("border", "none")))
       },
       fasteval=function(X,U,EqualPairs,pairpot,potpars,correction,
                         ..., halfway=FALSE) {
         # fast evaluator for Geyer interaction
         if(!all(correction %in% c("border", "none")))
           return(NULL)
         if(spatstat.options("fasteval") == "test")
           message("Using fast eval for Geyer")
         r  <- potpars$r
         hc <- potpars$sat
         # first determine saturated pair counts
         counts <- strausscounts(U, X, r, EqualPairs) 
         satcounts <- pmin(sat, counts)
         satcounts <- matrix(satcounts, ncol=1)
         if(halfway) {
           # trapdoor used by suffstat() 
           return(satcounts)
         } else if(sat == Inf) {
           # no saturation: fast code
           return(2 * satcounts)
         }
         # extract counts for data points
         Uindex <- EqualPairs[,2]
         Xindex <- EqualPairs[,1]
         Xcounts <- integer(npoints(X))
         Xcounts[Xindex] <- counts[Uindex]
         # evaluate change in saturated counts of other data points
         change <- geyercounts(U, X, r, sat, Xcounts, EqualPairs)
         answer <- satcounts + change
         return(matrix(answer, ncol=1))
       }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}

geyercounts <- function(U, X, r, sat, Xcounts, EqualPairs) {
  # evaluate effect of adding dummy point or deleting data point
  # on saturated counts of other data points
  stopifnot(is.numeric(r))
  stopifnot(is.numeric(sat))
  # for C calls we need finite numbers
  stopifnot(is.finite(r))
  stopifnot(is.finite(sat))
  # sort in increasing order of x coordinate
  oX <- order(X$x)
  oU <- order(U$x)
  Xsort <- X[oX]
  Usort <- U[oU]
  nX <- npoints(X)
  nU <- npoints(U)
  Xcountsort <- Xcounts[oX]
  # inverse: data point i has sorted position i' = rankX[i]
  rankX <- integer(nX)
  rankX[oX] <- seq_len(nX)
  rankU <- integer(nU)
  rankU[oU] <- seq_len(nU)
  # map from quadrature points to data points
  Uindex <- EqualPairs[,2]
  Xindex <- EqualPairs[,1]
  Xsortindex <- rankX[Xindex]
  Usortindex <- rankU[Uindex]
  Cmap <- rep(-1, nU)
  Cmap[Usortindex] <- Xsortindex - 1
  # call C routine
  zz <- .C("Egeyer",
           nnquad = as.integer(nU),
           xquad  = as.double(Usort$x),
           yquad  = as.double(Usort$y),
           quadtodata = as.integer(Cmap),
           nndata = as.integer(nX),
           xdata  = as.double(Xsort$x),
           ydata  = as.double(Xsort$y),
           tdata  = as.integer(Xcountsort),
           rrmax  = as.double(r),
           ssat   = as.double(sat),
           result = as.double(numeric(nU)),
           PACKAGE="spatstat")
  result <- zz$result[rankU]
  return(result)
}

