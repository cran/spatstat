#
# split.R
#
# $Revision: 1.4 $ $Date: 2006/02/25 02:17:22 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = x$marks, drop=FALSE, un=NULL, ...) {
  verifyclass(x, "ppp")
  if(is.null(un))
     un <- missing(f)
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
  lev <- levels(f)
  if(drop) lev <- lev[table(f) > 0]
  out <- list()
  for(l in lev) 
    out[[paste(l)]] <- x[f == l]
  
  if(un)
     out <- lapply(out, unmark)
  class(out) <- c("splitppp", class(out))
  return(out)
}

"split<-.ppp" <- function(x, f=x$marks, drop=FALSE, un=missing(f), 
                          ..., value) {
  verifyclass(x, "ppp")
  stopifnot(is.list(value))
  if(!all(unlist(lapply(value, is.ppp))))
    stop("Each entry of \`value\' must be a point pattern")

  ismark <- unlist(lapply(value, is.marked))
  if(any(ismark) && !all(ismark))
    stop("Some entries of \`value\' are marked, and others are unmarked")
  vmarked <- all(ismark)

  # evaluate `un' before assigning value of 'f'
  un <- un
  
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
  if(!drop) {
    lev <- levels(f)
    levtype <- "levels of f"
  } else {
    lev <- levels(f)[table(f) > 0]
    levtype <- "levels which f actually takes"
  }
  if(length(value) != length(lev))
      stop(paste("length of \`value\' should equal the number of",
                 levtype))

  # ensure value[[i]] is associated with lev[i]
  if(!is.null(names(value))) {
    if(!all(names(value) %in% paste(levels(f))))
      stop("names of \`value\' should be levels of f")
    value <- value[lev]
  }
  names(value) <- NULL
  
  # restore the marks, if they were discarded
  if(un && is.marked(x)) {
    if(vmarked)
      warning("\`value\' contains marked point patterns:\
 this is inconsistent with un=TRUE; marks ignored.")
    for(i in seq(value)) 
      value[[i]] <- value[[i]] %mark% factor(lev[i], levels=levels(f))
  }

  out <- superimpose(value)
  return(out)
}
