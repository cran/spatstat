#
# split.ppp.R
#
# $Revision: 1.9 $ $Date: 2008/07/22 21:15:54 $
#
# split.ppp and "split<-.ppp"
#
#########################################

split.ppp <- function(x, f = marks(x), drop=FALSE, un=NULL, ...) {
  verifyclass(x, "ppp")
  if(is.null(un))
    un <- missing(f)

  if(missing(f)) {
    # f defaults to marks of x
    if(!is.multitype(x)) 
      stop("f is missing and there is no sensible default")
    f <- fsplit <- marks(x)
  } else{
    # f was given
    fsplit <- f
    # if f is a tessellation or image, determine the grouping
    if(inherits(f, "im"))
      fsplit <- f <- tess(image=f)
    if(inherits(f, "tess")) {
      f <- marks(cut(x, f))
    } else {
      if(!is.factor(f))
        stop("f must be a factor")
      if(length(f) != x$n)
        stop("length(f) must equal the number of points in x")
    }
  }

  lev <- levels(f)
  if(drop) {
    lev <- lev[table(f) > 0]
    fsplit <- fsplit[f %in% lev]
  }

  # split the data
  out <- list()
  for(l in lev) 
    out[[paste(l)]] <- x[f == l]
  
  if(un)
     out <- lapply(out, unmark)
  if(inherits(fsplit, "tess")) {
    til <- tiles(fsplit)
    for(i in seq(along=out))
      out[[i]]$window <- til[[i]]
  }
  class(out) <- c("splitppp", class(out))
  attr(out, "fsplit") <- fsplit
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
    fsplit <- f
    if(inherits(f, "tess"))
      f <- marks(cut(x, f))
    if(!is.factor(f))
      stop("f must be a factor")
    if(length(f) != x$n)
      stop("length(f) must equal the number of points in x")
  } else {
    if(is.multitype(x))
      f <- fsplit <- marks(x)
    else
      stop("f is missing and there is no sensible default")
  }
  lev <- levels(f)
  if(!drop) 
    levtype <- "levels of f"
  else {
    levtype <- "levels which f actually takes"
    lev <- lev[table(f) > 0]
    fsplit <- fsplit[f %in% lev]
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


print.splitppp <- function(x, ...) {
  f <- attr(x, "fsplit")
  cat(paste("Point pattern split by",
            if(inherits(f, "tess")) "tessellation" else "factor",
            "\n"))
  nam <- names(x)
  for(i in seq(length(x))) {
    cat(paste("\n", nam[i], ":\n", sep=""))
    print(x[[i]])
  }
  return(invisible(NULL))
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

"[.splitppp" <- function(x, ...) {
  f <- attr(x, "fsplit")
  # invoke list method on x
  class(x) <- "list"
  y <- x[...]
  # then make it a 'splitppp' object too
  class(y) <- c("splitppp", class(y))
  attr(y, "fsplit") <- f[...]
  y
}

"[<-.splitppp" <- function(x, ..., value) {
  f <- attr(x, "fsplit")
  # invoke list method
  class(x) <- "list"
  x[...] <- value
  # then make it a 'splitppp' object too
  class(x) <- c("splitppp", class(x))
  attr(x, "fsplit") <- f
  x
}
  
