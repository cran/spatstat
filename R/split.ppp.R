#
# split.R
#
# $Revision: 1.3 $ $Date: 2005/08/10 06:55:02 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = x$marks) {
  verifyclass(x, "ppp")
  if(!missing(f)) {
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(is.marked(x) && is.factor(x$marks)) 
      f <- x$marks
    else
      stop("f is missing and there is no sensible default")
  }
  out <- list()
  for(l in levels(f)) 
    out[[paste(l)]] <- x[f == l]
  
  class(out) <- c("splitppp", class(out))
  return(out)
}

"split<-.ppp" <- function(x, f=x$marks, value) {
  verifyclass(x, "ppp")
  fimplicit <- FALSE
  if(!missing(f)) {
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(!is.marked(x) || !is.factor(x$marks))
      stop("f is missing and there is no sensible default")
    f <- x$marks
    fimplicit <- TRUE
  }
  if(!is.list(value) || length(value) != length(levels(f)))
    stop("value must be a list, with an entry for each level of f")
  if(!all(unlist(lapply(value, is.ppp))))
    stop("Each entry of \`value\' must be a point pattern")

  if(is.null(names(value))) {
    if(fimplicit)
      names(value) <- paste(levels(f)) # these become the marks
  }else {
    if(!all(names(value) %in% paste(levels(f))))
      stop("names of \`value\' should be levels of f")
  }

  out <- superimpose(value)

  return(out)
}
