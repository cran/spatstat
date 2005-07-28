#
#
#    multipair.util.R
#
#    $Revision: 1.7 $	$Date: 2005/07/27 06:32:16 $
#
#    Utilities for multitype pairwise interactions
#	
# -------------------------------------------------------------------
#	


MultiPair.checkmatrix <-
  function(mat, n, name) {
    if(!is.matrix(mat))
      stop(paste(name, "must be a matrix"))
    if(any(dim(mat) != rep(n,2)))
      stop(paste(name, "must be a square matrix,",
                 "of size", n, "x", n))
    isna <- is.na(mat)
    if(any(mat[!isna] <= 0))
      stop(paste("Entries of", name,
                 "must be positive numbers or NA"))
    if(any(isna != t(isna)) ||
       any(mat[!isna] != t(mat)[!isna]))
      stop(paste(name, "must be a symmetric matrix"))
  }

