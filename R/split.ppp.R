#
# split.R
#
# $Revision: 1.1 $ $Date: 2004/08/26 05:12:58 $
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
  if(!is.list(value) || length(value) != length(levels(f)))
    stop("value must be a list, with an entry for each level of f")
  if(!all(unlist(lapply(value, is.ppp))))
    stop("Each entry of \`value\' must be a point pattern")

  if(is.null(names(value)))
    names(value) <- paste(levels(f))
  
  out <- x
  for(l in levels(f))
    out[f == l] <- value[[paste(l)]]

  return(out)
}
