#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.2 $   $Date: 2006/10/10 04:22:48 $
#
#
rlabel <- function(X, labels=marks(X), permute=TRUE) {
  verifyclass(X, "ppp")
  if(is.null(labels))
    stop("labels not given and marks not present")
  Y <- X %mark% sample(labels, X$n, replace=!permute)
  return(Y)
}

