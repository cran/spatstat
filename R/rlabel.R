#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.4 $   $Date: 2010/03/08 12:46:35 $
#
#
rlabel <- function(X, labels=marks(X), permute=TRUE) {
  verifyclass(X, "ppp")
  if(is.null(labels))
    stop("labels not given and marks not present")
  npoints <- X$n
  if(is.vector(labels) || is.factor(labels)) {
    nlabels <- length(labels)
    if(permute && (nlabels != npoints))
      stop("length of labels vector does not match number of points")
    Y <- X %mark% sample(labels, npoints, replace=!permute)
  } else if(is.data.frame(labels)) {
    nlabels <- nrow(labels)
    if(permute && (nlabels != npoints))
      stop("number of rows of data frame does not match number of points")      
    Y <- X %mark% labels[sample(1:nlabels, npoints, replace=!permute), ]
  } else stop("Format of labels argument is not understood")
  return(Y)
}

