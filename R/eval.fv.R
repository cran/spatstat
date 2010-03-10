#
#     eval.fv.R
#
#
#        eval.fv()             Evaluate expressions involving fv objects
#
#        compatible.fv()       Check whether two fv objects are compatible
#
#     $Revision: 1.10 $     $Date: 2010/01/08 23:28:06 $
#

eval.fv <- function(expr, envir) {
  # convert syntactic expression to 'expression' object
  e <- as.expression(substitute(expr))
  # convert syntactic expression to call
  elang <- substitute(expr)
  # find names of all variables in the expression
  varnames <- all.vars(e)
  if(length(varnames) == 0)
    stop("No variables in this expression")
  # get the actual variables
  if(missing(envir))
    envir <- sys.parent()
  vars <- lapply(as.list(varnames), function(x, ee) get(x, envir=ee), ee=envir)
  names(vars) <- varnames
  # find out which ones are fv objects
  fvs <- unlist(lapply(vars, is.fv))
  if(!any(fvs))
    stop("No fv objects in this expression")
  funs <- vars[fvs]
  nfuns <- length(funs)
  # test whether the fv objects are compatible
  if(nfuns > 1) {
    # test compatibility
    for(i in 2:nfuns)
      if(!compatible.fv(funs[[1]], funs[[i]]))
        stop(paste("Objects", names(funs)[1], "and", names(funs)[i],
                   "are incompatible"))
  }
  # copy first object as template
  result <- funs[[1]]
  # determine which function estimates are supplied
  argname <- fvnames(result, ".x")
  nam <- names(result)
  ynames <- nam[nam != argname]
  # for each function estimate, evaluate expression
  for(yn in ynames) {
    # extract corresponding estimates from each fv object
    funvalues <- lapply(funs, function(x, n) { x[[n]] }, n=yn)
    # insert into list of argument values
    vars[fvs] <- funvalues
    # evaluate
    result[[yn]] <- eval(e, vars)
  }
  # determine y axis label for the result
  ylabs <- lapply(funs, function(x) { attr(x, "ylab") })
  attr(result, "ylab") <-
    eval(substitute(substitute(e, ylabs), list(e=elang)))
  attr(result, "fname") <- "eval"
  return(result)
}

compatible.fv <- function(A, B) {
  verifyclass(A, "fv")
  verifyclass(B, "fv")
  namesmatch <-
    identical(all.equal(names(A),names(B)), TRUE) &&
    (fvnames(A, ".x") == fvnames(B, ".x")) &&
    (fvnames(A, ".y") == fvnames(B, ".y"))
  rA <- with(A, .x)
  rB <- with(B, .x)
  approx.equal <- function(x, y) { max(abs(x-y)) <= .Machine$double.eps }
  rmatch <- (length(rA) == length(rB)) && approx.equal(rA, rB)
  return(namesmatch && rmatch)
}

