#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.1 $   $Date: 2005/06/08 01:33:39 $
#
#
rlabel <- function(X, labels=X$marks, permute=TRUE) {
  verifyclass(X, "ppp")
  if(is.null(labels))
    stop("labels not given and marks not present")
  Y <- X %mark% sample(labels, X$n, replace=!permute)
  return(Y)
}

