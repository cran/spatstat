#
#     eval.im.R
#
#        eval.im()             Evaluate expressions involving images
#
#        compatible.im()       Check whether two images are compatible
#
#     $Revision: 1.16 $     $Date: 2010/07/16 04:18:05 $
#

eval.im <- function(expr, envir) {
  e <- as.expression(substitute(expr))
  # get names of all variables in the expression
  varnames <- all.vars(e)
  allnames <- all.names(e, unique=TRUE)
  funnames <- allnames[!(allnames %in% varnames)]
  if(length(varnames) == 0)
    stop("No variables in this expression")
  # get the values of the variables
  if(missing(envir))
    envir <- sys.parent()
  vars <- lapply(as.list(varnames), function(x, e) get(x, envir=e), e=envir)
  names(vars) <- varnames
  funs <- lapply(as.list(funnames), function(x, e) get(x, envir=e), e=envir)
  names(funs) <- funnames
  # find out which variables are images
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
    v <- as.matrix(x)
    dim(v) <- NULL
    return(v)
  }
  imagevalues <- lapply(images, getvalues)
  template <- images[[1]]
  # This bit has been repaired:
  vars[ims] <- imagevalues
  v <- eval(e, append(vars, funs))
  #
  # reshape, etc
  result <- im(v, template$xcol, template$yrow, 
               unitname=unitname(template))
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
  uok <- identical(all.equal(unitname(A), unitname(B)), TRUE)
  return(xok && yok && uok)
}

