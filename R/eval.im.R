#
#     eval.im.R
#
#        eval.im()             Evaluate expressions involving images
#
#        compatible.im()       Check whether two images are compatible
#
#     $Revision: 1.3 $     $Date: 2005/07/28 06:11:44 $
#

eval.im <- function(expr) {
  e <- as.expression(substitute(expr))
  # get names of all variables in the expression
  varnames <- all.vars(e)
  if(length(varnames) == 0)
    stop("No variables in this expression")
  # get the values of the variables
  pe <- sys.parent()
  vars <- lapply(as.list(varnames), function(x, e) get(x, envir=e), e=pe)
  names(vars) <- varnames
  # find out which are images
  ims <- unlist(lapply(vars, is.im))
  if(!any(ims))
    stop("No images in this expression")
  images <- vars[ims]
  nimages <- length(images)
  # test the images are compatible
  if(nimages > 1) {
    # test compatibility
    for(i in 2:nimages)
      if(!compatible.im(images[[1]], images[[i]]))
        stop(paste("Images", names(images)[1], "and", names(images)[i],
                   "are incompatible"))
  }
  # replace each image by its matrix of pixel values, and evaluate 
  imagevalues <- lapply(images, function(x) { x$v })
  result <- images[[1]]
  v <- eval(e, imagevalues)
  # sanity check
  if(!is.matrix(v))
    stop("Evaluating the expression did not yield a matrix")
  if(any(dim(v) != dim(result$v)))
     stop("Expression yields a matrix with inappropriate dimensions")
  result$v <- v
  return(result)
}
  
compatible.im <- function(A, B) {
  verifyclass(A, "im")
  verifyclass(B, "im")
  agree <- function(x, y) { max(abs(x-y)) <= .Machine$double.eps }
  return((all(A$dim == B$dim)) &&
         agree(A$xrange, B$xrange) &&
         agree(A$yrange, B$yrange) &&
         agree(A$xstep, B$xstep) &&
         agree(A$ystep, B$ystep) &&
         agree(A$xcol, B$xcol) &&
         agree(A$yrow, B$yrow))
}

