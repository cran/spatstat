#
#
#    multistrauss.S
#
#    $Revision: 1.7 $	$Date: 2001/11/05 08:20:43 $
#
#    The multitype Strauss process
#
#    MultiStrauss()    create an instance of the multitype Strauss process
#                 [an object of class 'interact']
#	
# -------------------------------------------------------------------
#	

MultiStrauss <- function(types, radii) {
  if(length(types) == 1)
    stop("The \`types\' argument should be a vector of all possible types")
  if(!is.factor(types)) {
    # warning("\`types\' converted to a factor")
    types <- factor(types)
  }
  ct <- levels(types)
  dimnames(radii) <- list(ct, ct)
  out <- 
  list(
         name     = "Multitype Strauss process",
         family    = multipair.family,
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

     # list all UNORDERED pairs of types to be checked
     # (the interaction must be symmetric in type, and scored as such)
     uptri <- (row(r) <= col(r)) & !is.na(r)
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
     # assemble the relevant interaction distance for each pair of points
     rxu <- r[ tx, tu ]
     # apply relevant threshold to each pair of points
     str <- (d <= rxu)
     # create logical array for result
     z <- array(F, dim=c(dim(d), npairs),
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
         init     = function(self) {
                      r <- self$par$radii
                      nt <- length(self$par$types)
                      MultiPair.checkmatrix(r, nt, "\`radii\'")
                    },
         update = NULL,  # default OK
         print = function(self) {
           print.isf(self$family)
           cat(paste("Interaction:\t", self$name, "\n"))
           cat(paste(length(self$par$types), "types of points\n"))
           cat("Possible types: \n")
           print(self$par$types)
           cat("Interaction radii:\n")
           print(self$par$radii)
           invisible()
         }
  )
  class(out) <- "interact"
  out$init(out)
  return(out)
}
