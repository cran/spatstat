#
#    util.S    miscellaneous utilities
#
#    $Revision: 1.4 $    $Date: 2002/04/07 11:11:54 $
#
#  (a) for matrices only:
#
#    matrowany(X) is equivalent to apply(X, 1, any)
#    matrowall(X) "   "  " "  "  " apply(X, 1, all)
#    matcolany(X) "   "  " "  "  " apply(X, 2, any)
#    matcolall(X) "   "  " "  "  " apply(X, 2, all)
#
#  (b) for 3D arrays only:
#    apply23sum(X)  "  "   "  " apply(X, c(2,3), sum)
#
#  (c) weighted histogram
#    whist()
#
matrowsum <- function(x) {
  x %*% rep(1, ncol(x))
}

matcolsum <- function(x) {
  rep(1, nrow(x)) %*% x
}
  
matrowany <- function(x) {
  (matrowsum(x) > 0)
}

matrowall <- function(x) {
  (matrowsum(x) == ncol(x))
}

matcolany <- function(x) {
  (matcolsum(x) > 0)
}

matcolall <- function(x) {
  (matcolsum(x) == nrow(x))
}

########
    # hm, this is SLOWER

apply23sum <- function(x) {
  dimx <- dim(x)
  if(length(dimx) != 3)
    stop("x is not a 3D array")
  result <- array(0, dimx[-1])

  nz <- dimx[3]
  for(k in 1:nz) {
    result[,k] <- matcolsum(x[,,k])
  }
  result
}
    
#######################
#
#    whist      weighted histogram
#

whist <-
  function(x, breaks, weights) {
    if(missing(weights)) 
      h <- hist(x, breaks=breaks, plot=FALSE,probability=FALSE)$counts
    else {
      # Thanks to Peter Dalgaard
      cell <- cut(x, breaks, include.lowest=TRUE)
      h <- tapply(weights, cell, sum)
      h[is.na(h)] <- 0
    }
    return(h)
}
