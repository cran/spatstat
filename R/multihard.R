#
#
#    multihard.R
#
#    $Revision: 1.2 $	$Date: 2011/05/17 13:21:11 $
#
#    The Hard core process
#
#    Hardcore()     create an instance of the Hard Core process
#                      [an object of class 'interact']
#	
#
# -------------------------------------------------------------------
#	

MultiHard <- function(types=NA, hradii) {
  nt <- length(types)
  if(nt > 1) {
    if(is.factor(types)) {
      types <- levels(types)
    } else {
      types <- levels(factor(types, levels=types))
    }
    dimnames(hradii) <- list(types, types)
  } else if((nt == 0) || !is.na(types))
    stop(paste("The", sQuote("types"),"argument should",
                           "either be a vector of all possible types or \"NA\".\n"))
  out <- 
  list(
         name     = "Multitype Hardcore process",
         creator  = "MultiHard",
         family    = pairwise.family,
         pot      = function(d, tx, tu, par) {
     # arguments:
     # d[i,j] distance between points X[i] and U[j]
     # tx[i]  type (mark) of point X[i]
     # tu[i]  type (mark) of point U[j]
     #
     # get matrices of interaction radii
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
     uptri <- (row(h) <= col(h)) & (!is.na(h))
     mark1 <- (lx[row(h)])[uptri]
     mark2 <- (lx[col(h)])[uptri]
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
     # apply the relevant hard core distance to each pair of points
     hxu <- h[ tx, tu ]
     forbid <- (d < hxu)
     forbid[is.na(forbid)] <- FALSE
     # form the potential 
     value <- ifelse(forbid, -Inf, 0)
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
     attr(z, "IsOffset") <- TRUE
     return(z)
     },
     #### end of 'pot' function ####
     #       
         par      = list(types=types, hradii = hradii),
         parnames = c("possible types", "hardcore distances"),
         selfstart = function(X, self) {
		if(length(self$par$types) > 1) return(self)
                types <- levels(marks(X))
                MultiHard(types=types,hradii=self$par$hradii)
	 },
         init     = function(self) {
                      nt <- length(self$par$types)
                      chk <- (nt != 1)# || is.na(self$par$types)
                      if(chk) {
                        h <- self$par$hradii
                        MultiPair.checkmatrix(h, nt, sQuote("hradii"))
                      }
                    },
         update = NULL,  # default OK
         print = function(self) {
           print.isf(self$family)
           cat(paste("Interaction:\t", self$name, "\n"))
           cat(paste(length(self$par$types), "types of points\n"))
           cat("Possible types: \n")
           print(self$par$types)
           cat("Hardcore radii:\n")
           print(self$par$hradii)
           invisible()
         },
        interpret = function(coeffs, self) {
          # there are no regular parameters (woo-hoo!)
          return(NULL)
        },
        valid = function(coeffs, self) {
          return(TRUE)
        },
        project = function(coeffs, self) {
          return(coeffs)
        },
         irange = function(self, coeffs=NA, epsilon=0, ...) {
           h <- self$par$hradii
           return(max(0, h, na.rm=TRUE))
         },
       version=versionstring.spatstat()
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
