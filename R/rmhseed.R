#
#  rmhseed.R
#
#  $Revision: 1.1 $   $Date: 2005/03/10 02:51:45 $
#

rmhseed <- function(iseed = NULL) {
  if(is.null(iseed))
    iseed <- sample(1:1000000,3)
  else if(!is.numeric(iseed) || length(iseed) != 3)
    stop("iseed should be a vector of 3 integers")
  iseed.save <- iseed

# The R in-house random number/sampling system gets used in rmh.default.
# Therefore we build an argument for set.seed() from iseed
# which determines all the randomness, including ``proposal points''
# and possibly the starting state.
# Thus supplying the same ``iseed'' in a repeat simulation will
# give ***exactly*** the same result.

  rrr <- .Fortran(
                  "arand",
                  ix=as.integer(iseed[1]),
                  iy=as.integer(iseed[2]),
                  iz=as.integer(iseed[3]),
                  rand=double(1),
                  PACKAGE="spatstat"
                  )
  build.seed <- round(rrr$rand*1e6)
  iseed <- unlist(rrr[c("ix","iy","iz")])

# return seed object  
  out <- list(iseed=iseed,
              iseed.save=iseed.save,
              build.seed=build.seed)
  class(out) <- c("rmhseed", class(out))
  return(out)
}


print.rmhseed <- function(x, ...) {
   cat("Random number seeds:\n")
   cat(paste("iseed (original) = (",
             paste(x$iseed.save, collapse=", "),
             ")\n", sep=""))
   cat(paste("iseed (computed) = (",
             paste(x$iseed, collapse=", "),
             ")\n", sep=""))
   cat(paste("build seed = ", x$build.seed, "\n"))
}
