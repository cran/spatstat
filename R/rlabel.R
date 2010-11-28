#
#   rlabel.R
#
#   random (re)labelling
#
#   $Revision: 1.5 $   $Date: 2010/11/21 04:24:29 $
#
#
rlabel <- function(X, labels=marks(X), permute=TRUE) {
  verifyclass(X, "ppp")
  if(is.null(labels))
    stop("labels not given and marks not present")
  npts <- npoints(X)
  if(is.vector(labels) || is.factor(labels)) {
    nlabels <- length(labels)
    if(permute && (nlabels != npts))
      stop("length of labels vector does not match number of points")
    Y <- X %mark% sample(labels, npts, replace=!permute)
  } else if(is.data.frame(labels)) {
    nlabels <- nrow(labels)
    if(permute && (nlabels != npts))
      stop("number of rows of data frame does not match number of points")      
    Y <- X %mark% labels[sample(1:nlabels, npts, replace=!permute), ]
  } else stop("Format of labels argument is not understood")
  return(Y)
}

