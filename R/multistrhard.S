#
#
#    multistrhard.S
#
#    $Revision: 2.16 $	$Date: 2007/01/11 03:36:02 $
#
#    The multitype Strauss/hardcore process
#
#    MultiStraussHard()
#                 create an instance of the multitype Strauss/ harcore
#                 point process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

MultiStraussHard <- function(types, iradii, hradii) {
  if(length(types) == 1)
    stop(paste("The", sQuote("types"),
               "argument should be a vector of all possible types"))
  if(is.factor(types)) {
    types <- levels(types)
  } else {
    types <- levels(factor(types, levels=types))
  }
  dimnames(iradii) <- dimnames(hradii) <- list(types, types)
  out <- 
  list(
         name     = "Multitype Strauss Hardcore process",
         creator  = "MultiStraussHard",
         family    = pairwise.family,
         pot      = function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[i]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
     r <- par$iradii
     h <- par$hradii

     # get possible marks and validate
     if(!is.factor(tx) || !is.factor(tu))
	stop("marks of data and dummy points must be factor variables")
     lx <- levels(tx)
     lu <- levels(tu)
     if(length(lx) != length(lu) || any(lx != lu))
	stop("marks of data and dummy points do not have same possible levels")

     if(!identical(lx, par$types))
        stop("data and model do not have the same possible levels of marks")
     if(!identical(lu, par$types))
        stop("dummy points and model do not have the same possible levels of marks")
                   
     # list all UNORDERED pairs of types to be checked
     # (the interaction must be symmetric in type, and scored as such)
     uptri <- (row(r) <= col(r)) & (!is.na(r) | !is.na(h))
     mark1 <- (lx[row(r)])[uptri]
     mark2 <- (lx[col(r)])[uptri]
     vname <- apply(cbind(mark1,mark2), 1, paste, collapse="x")
     vname <- paste("mark", vname, sep="")
     vname <- make.names(vname)  # converts illegal characters
     npairs <- length(vname)
     # list all ORDERED pairs of types to be checked
     # (to save writing the same code twice)
     different <- mark1 != mark2
     mark1o <- c(mark1, mark2[different])
     mark2o <- c(mark2, mark1[different])
     nordpairs <- length(mark1o)
     # unordered pair corresponding to each ordered pair
     ucode <- c(1:npairs, (1:npairs)[different])
     #
     # go....
     # apply the relevant interaction distance to each pair of points
     rxu <- r[ tx, tu ]
     str <- (d < rxu)
     str[is.na(str)] <- FALSE
     # and the relevant hard core distance
     hxu <- h[ tx, tu ]
     forbid <- (d < hxu)
     forbid[is.na(forbid)] <- FALSE
     # form the potential 
     value <- ifelse(forbid, -Inf, str)
     # create numeric array for result
     z <- array(0, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # assign value[i,j] -> z[i,j,k] where k is relevant interaction code
     for(i in 1:nordpairs) {
       # data points with mark m1
       Xsub <- (tx == mark1o[i])
       # quadrature points with mark m2
       Qsub <- (tu == mark2o[i])
       # assign
       z[Xsub, Qsub, ucode[i]] <- value[Xsub, Qsub]
     }     
     return(z)
     },
     #### end of 'pot' function ####
     #       
         par      = list(types=types, iradii = iradii, hradii = hradii),
         parnames = c("possible types", "interaction distances", "hardcore distances"),
         init     = function(self) {
                      r <- self$par$iradii
                      h <- self$par$hradii
                      nt <- length(self$par$types)

                      MultiPair.checkmatrix(r, nt, sQuote("iradii"))
                      MultiPair.checkmatrix(h, nt, sQuote("hradii"))

                      ina <- is.na(iradii)
                      hna <- is.na(hradii)
                      if(all(ina))
                        stop(paste("All entries of", sQuote("iradii"),
                                   "are NA"))
                      both <- !ina & !hna
                      if(any(iradii[both] <= hradii[both]))
                        stop("iradii must be larger than hradii")
                    },
         update = NULL,  # default OK
         print = function(self) {
           print.isf(self$family)
           cat(paste("Interaction:\t", self$name, "\n"))
           cat(paste(length(self$par$types), "types of points\n"))
           cat("Possible types: \n")
           print(self$par$types)
           cat("Interaction radii:\n")
           print(self$par$iradii)
           cat("Hardcore radii:\n")
           print(self$par$hradii)
           invisible()
         },
        interpret = function(coeffs, self) {
          # get possible types
          typ <- self$par$types
          ntypes <- length(typ)
          # get matrices of interaction radii
          r <- self$par$iradii
          h <- self$par$hradii
          # list all relevant unordered pairs of types
          uptri <- (row(r) <= col(r)) & (!is.na(r) | !is.na(h))
          index1 <- (row(r))[uptri]
          index2 <- (col(r))[uptri]
          npairs <- length(index1)
          # extract canonical parameters; shape them into a matrix
          gammas <- matrix(, ntypes, ntypes)
          dimnames(gammas) <- list(typ, typ)
          expcoef <- exp(coeffs)
          gammas[ cbind(index1, index2) ] <- expcoef
          gammas[ cbind(index2, index1) ] <- expcoef
          #
          return(list(param=list(gammas=gammas),
                      inames="interaction parameters gamma_ij",
                      printable=round(gammas,4)))
        },
        valid = function(coeffs, self) {
           # interaction radii r[i,j]
           iradii <- self$par$iradii
           # hard core radii r[i,j]
           hradii <- self$par$hradii
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # Check that we managed to estimate all required parameters
           required <- !is.na(iradii)
           if(!all(is.finite(gamma[required])))
             return(FALSE)
           # Check that the model is integrable
           # inactive hard cores ...
           ihc <- (is.na(hradii) | hradii == 0)
           # .. must have gamma <= 1
           return(all(gamma[required & ihc] <= 1))
         },
         project = function(coeffs, self) {
           # interaction radii r[i,j]
           iradii <- self$par$iradii
           # hard core radii r[i,j]
           hradii <- self$par$hradii
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # remove NA's
           gamma[is.na(gamma)] <- 1
           # inactive hard cores?
           ihc <- (is.na(hradii) | hradii == 0)
           if(any(ihc & (gamma > 1))) 
               gamma[ihc] <- pmin(gamma[ihc], 1)
           # now put them back..
           # get matrices of interaction radii
           r <- self$par$iradii
           h <- self$par$hradii
           # list all unordered pairs of types
           uptri <- (row(r) <= col(r)) & (!is.na(r) | !is.na(h))
           # reassign 
           coeffs[] <- log(gamma[uptri])
           return(coeffs)
        },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$iradii
           h <- self$par$hradii
           ractive <- !is.na(r)
           hactive <- !is.na(h)
           if(any(!is.na(coeffs))) {
             gamma <- (self$interpret)(coeffs, self)$param$gammas
             gamma[is.na(gamma)] <- 1
             ractive <- ractive & (abs(log(gamma)) > epsilon)
           }
           if(!any(c(ractive,hactive)))
             return(0)
           else
             return(max(c(r[ractive],h[hactive])))
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
