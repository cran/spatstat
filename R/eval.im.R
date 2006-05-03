#
#     eval.im.R
#
#        eval.im()             Evaluate expressions involving images
#
#        compatible.im()       Check whether two images are compatible
#
#     $Revision: 1.10 $     $Date: 2006/05/01 11:19:10 $
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
  getvalues <- function(x) {
    v <- as.vector(x$v)
    if(x$type != "factor") return(v)
    else return(factor(v, levels=seq(x$lev), labels=lev))
  }
  imagevalues <- lapply(images, getvalues)
  template <- images[[1]]
  # This bit has been repaired:
  vars[ims] <- imagevalues
  v <- eval(e, vars)
  #
  # reshape, etc
  result <- im(v, template$xcol, template$yrow, levels(v))
  return(result)
}
  
compatible.im <- function(A, B, tol=1e-6) {
  verifyclass(A, "im")
  verifyclass(B, "im")
  xdiscrep <- max(abs(A$xrange - B$xrange),
                 abs(A$xstep - B$xstep),
                 abs(A$xcol - B$xcol))
  ydiscrep <- max(abs(A$yrange - B$yrange),
                 abs(A$ystep - B$ystep),
                 abs(A$yrow - B$yrow))
  xok <- (xdiscrep < tol * min(A$xstep, B$xstep))
  yok <- (ydiscrep < tol * min(A$ystep, B$ystep))
  return(xok && yok)
}

