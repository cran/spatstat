#
#
#    multistrauss.S
#
#    $Revision: 2.15 $	$Date: 2009/10/22 20:41:37 $
#
#    The multitype Strauss process
#
#    MultiStrauss()    create an instance of the multitype Strauss process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

MultiStrauss <- function(types=NULL, radii) {
  if(!is.null(types)) {
    if(length(types) == 0)
      stop(paste("The", sQuote("types"),"argument should be",
                 "either NULL or a vector of all possible types"))
    if(any(is.na(types)))
      stop("NA's not allowed in types")
    if(is.factor(types)) {
      types <- levels(types)
    } else {
      types <- levels(factor(types, levels=types))
    }
    dimnames(radii) <- list(types, types)
  } 
  out <- 
  list(
         name     = "Multitype Strauss process",
         creator  = "MultiStrauss",
         family    = pairwise.family,
         pot      = function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[j]  type (mark) of point U[j]
     #
     # get matrix of interaction radii r[ , ]
     r <- par$radii
     #
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
     uptri <- (row(r) <= col(r)) & !is.na(r)
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
     # assemble the relevant interaction distance for each pair of points
     rxu <- r[ tx, tu ]
     # apply relevant threshold to each pair of points
     str <- (d <= rxu)
     # create logical array for result
     z <- array(FALSE, dim=c(dim(d), npairs),
                dimnames=list(character(0), character(0), vname))
     # assign str[i,j] -> z[i,j,k] where k is relevant interaction code
     for(i in 1:nordpairs) {
       # data points with mark m1
       Xsub <- (tx == mark1o[i])
       # quadrature points with mark m2
       Qsub <- (tu == mark2o[i])
       # assign
       z[Xsub, Qsub, ucode[i]] <- str[Xsub, Qsub]
     }
     return(z)
   },
     #### end of 'pot' function ####
     #       
         par      = list(types=types, radii = radii),
         parnames = c("possible types", "interaction distances"),
         selfstart = function(X, self) {
		if(!is.null(self$par$types)) return(self)
                types <- levels(marks(X))
                MultiStrauss(types=types,radii=self$par$radii)
	 },
         init     = function(self) {
                      types <- self$par$types
                      if(!is.null(types)) {
                        h <- self$par$radii
                        nt <- length(types)
                        MultiPair.checkmatrix(h, nt, sQuote("radii"))
                      }
                    },
         update = NULL, # default OK
         print = function(self) {
           print.isf(self$family)
           cat(paste("Interaction:\t", self$name, "\n"))
           radii <- self$par$radii
           cat(paste(nrow(radii), "types of points\n"))
           types <- self$par$types
           if(!is.null(types)) {
             cat("Possible types: \n")
             print(types)
           } else cat("Possible types:\t not yet determined\n")
           cat("Interaction radii:\n")
           print(self$par$radii)
           invisible()
         },
        interpret = function(coeffs, self) {
          # get possible types
          typ <- self$par$types
          ntypes <- length(typ)
          # get matrix of Strauss interaction radii
          r <- self$par$radii
          # list all unordered pairs of types
          uptri <- (row(r) <= col(r)) & (!is.na(r))
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
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # interaction radii
           radii <- self$par$radii
           # parameters to estimate
           required <- !is.na(radii)
           gr <- gamma[required]
           return(all(is.finite(gr) & gr <= 1))
         },
         project  = function(coeffs, self) {
           # interaction parameters gamma[i,j]
           gamma <- (self$interpret)(coeffs, self)$param$gammas
           # remove NA's
           gamma[is.na(gamma)] <- 1
           # constrain them
           gamma <- matrix(pmin(gamma, 1),
                           nrow=nrow(gamma), ncol=ncol(gamma))
           # now put them back 
           # get matrix of Strauss interaction radii
           r <- self$par$radii
           # list all unordered pairs of types
           uptri <- (row(r) <= col(r)) & (!is.na(r))
           # reassign 
           coeffs[] <- log(gamma[uptri])
           return(coeffs)
         },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           r <- self$par$radii
           active <- !is.na(r)
           if(any(!is.na(coeffs))) {
             gamma <- (self$interpret)(coeffs, self)$param$gammas
             gamma[is.na(gamma)] <- 1
             active <- active & (abs(log(gamma)) > epsilon)
           }
           if(any(active)) return(max(r[active])) else return(0)
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
