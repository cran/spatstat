#
#
#    multistrhard.S
#
#    $Revision: 1.7 $	$Date: 2002/05/13 12:41:10 $
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
    stop("The \`types\' argument should be a vector of all possible types")
  if(!is.factor(types)) {
    # warning("\`types\' converted to a factor")
    types <- factor(types)
  }
  ct <- levels(types)
  dimnames(iradii) <- dimnames(hradii) <- list(ct, ct)
  out <- 
  list(
         name     = "Multitype Strauss Hardcore process",
         family    = multipair.family,
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

     # list all UNORDERED pairs of types to be checked
     # (the interaction must be symmetric in type, and scored as such)
     uptri <- (row(r) <= col(r)) & (!is.na(r) | !is.na(h))
     mark1 <- (lx[row(r)])[uptri]
     mark2 <- (lx[col(r)])[uptri]
     vname <- apply(cbind(mark1,mark2), 1, paste, collapse="x")
     vname <- paste("mark", vname, sep="")
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

                      MultiPair.checkmatrix(r, nt, "\`iradii\'")
                      MultiPair.checkmatrix(h, nt, "\`hradii\'")
                      
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
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
