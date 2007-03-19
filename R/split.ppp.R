#
# split.R
#
# $Revision: 1.8 $ $Date: 2007/03/16 08:33:10 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = marks(x), drop=FALSE, un=NULL, ...) {
  verifyclass(x, "ppp")
  if(is.null(un))
     un <- missing(f)
  if(!missing(f)) {
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(is.multitype(x))
      f <- marks(x)
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

"split<-.ppp" <- function(x, f=marks(x), drop=FALSE, un=missing(f), 
                          ..., value) {
  verifyclass(x, "ppp")
  stopifnot(is.list(value))
  if(!all(unlist(lapply(value, is.ppp))))
    stop(paste("Each entry of", sQuote("value"),
               "must be a point pattern"))

  ismark <- unlist(lapply(value, is.marked))
  if(any(ismark) && !all(ismark))
    stop(paste("Some entries of",
               sQuote("value"),
               "are marked, and others are unmarked"))
  vmarked <- all(ismark)

  # evaluate `un' before assigning value of 'f'
  un <- un
  
  if(!missing(f)) {
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(is.multitype(x))
      f <- marks(x)
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
      stop(paste("length of", sQuote("value"),
                 "should equal the number of",
                 levtype))

  # ensure value[[i]] is associated with lev[i]
  if(!is.null(names(value))) {
    if(!all(names(value) %in% paste(levels(f))))
      stop(paste("names of", sQuote("value"), "should be levels of f"))
    value <- value[lev]
  }
  names(value) <- NULL
  
  # restore the marks, if they were discarded
  if(un && is.marked(x)) {
    if(vmarked)
      warning(paste(sQuote("value"), "contains marked point patterns:",
                    "this is inconsistent with un=TRUE; marks ignored."))
    for(i in seq(value)) 
      value[[i]] <- value[[i]] %mark% factor(lev[i], levels=levels(f))
  }

  out <- superimpose(value, W=x$window)
  return(out)
}


summary.splitppp <- function(object, ...) {
  x <- lapply(object, summary, ...)
  class(x) <- "summary.splitppp"
  x
}

print.summary.splitppp <- function(x, ...) {
  class(x) <- "listof"
  print(x)
  invisible(NULL)
}

"[.splitppp" <- function(x, i) {
  # invoke list method
  class(x) <- "list"
  y <- x[i]
  # then make it a 'splitppp' object too
  class(y) <- c("splitppp", class(y))
  y
}

"[<-.splitppp" <- function(x, i, value) {
  # invoke list method
  class(x) <- "list"
  x[i] <- value
  # then make it a 'splitppp' object too
  class(x) <- c("splitppp", class(x))
  x
}
  
